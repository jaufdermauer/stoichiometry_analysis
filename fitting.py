from scipy.stats import gaussian_kde, binom
import numpy as np

class BGMatrixQuadrant:
    def __init__(self, r1, r2, c1, c2):
        self.r1 = r1
        self.r2 = r2
        self.c1 = c1
        self.c2 = c2

    def __str__(self):
        return "{0} {1} {2} {3}".format(self.r1, self.r2, self.c1, self.c2)

    def sum(self, matrix):
        return matrix[self.r1:self.r2, self.c1:self.c2].sum()

    def count(self, matrix):
        return float(matrix[self.r1:self.r2, self.c1:self.c2].size)

    def values(self, matrix):
        return matrix[self.r1:self.r2, self.c1:self.c2].flatten()

    def average(self, matrix):
        return self.sum(matrix) / self.count(matrix)

    def split(self):
        """
        :return:
        """
        q1 = BGMatrixQuadrant(self.r1, (self.r1 + self.r2) / 2, self.c1, (self.c1 + self.c2) / 2)
        q2 = BGMatrixQuadrant(self.r1, (self.r1 + self.r2) / 2, ((self.c1 + self.c2) / 2) + 1, self.c2)
        q3 = BGMatrixQuadrant(((self.r1 + self.r2) / 2) + 1, self.r2, self.c1, (self.c1 + self.c2) / 2)
        q4 = BGMatrixQuadrant(((self.r1 + self.r2) / 2) + 1, self.r2, ((self.c1 + self.c2) / 2) + 1, self.c2)

        return [q1, q2, q3, q4]


class BGMatrix:
    def __init__(self, matrix, min_width):
        self.matrix = matrix
        self.min_width = min_width

    def get_min_squares(self):
        width, height = self.matrix.shape
        q0 = BGMatrixQuadrant(0, width-1, 0, height-1)
        qs = q0.split()

        res = []
        for q in qs:
            res.append(self.recursive_search(q))

        return res

    def recursive_search(self, quadrant):
        if quadrant.r2 - quadrant.r1 <= self.min_width:
            return quadrant

        qs = quadrant.split()

        min_q = qs.pop(0)
        min_v = min_q.sum(self.matrix)

        for q in qs:
            v = q.sum(self.matrix)
            if v < min_v:
                min_v = v
                min_q = q

        return self.recursive_search(min_q)


def get_kde_values(values, num_samples=100, num_peaks=1):
    """
    Return the values correspoding to the Kernel Density Estimation
    :param values:
    :return:
    xs, ys
    """
    # Todo: Set a bandwidth according to num_peaks
    # Todo: Catch exception
    kde = gaussian_kde(values, 'silverman')
    ini = np.min(values)
    end = np.max(values)
    xs = np.linspace(ini, end, num_samples)
    return xs, kde.pdf(xs)


def get_peak_by_kde(values, num_samples=100, num_peaks=1):
    """
    Return the x of the maximum value of the kde
    Useful when estimating gaussians
    :param values:
    :param num_samples:
    :param num_peaks:
    :return:
    """
    xs, ys = get_kde_values(values, num_samples, num_peaks)
    return xs[np.argmax(ys)].flatten(), xs, ys


def get_bins_number(values):
    """
    Return the number of bins for values (Freedman-Diaconis rule)
    :param values:
    :return:
    """
    iqr = np.subtract(*np.percentile(values, [75, 25]))
    h = 2*iqr*np.power(len(values), -1.0/3.0)
    return int(np.ceil((np.max(values)-np.min(values))/h))


def get_outer_grid(x, y, width, outer_width):
    """
    Get the coordinates, as with meshgrid, but flattened, of a ring around a square ROI
    :param x: center of the ROI
    :param y: center of the ROI
    :param width: half width of the ROI, like a radius
    :param outer_width: thick of the outer ring
    :return:
    """
    Xin, Yin = np.meshgrid(np.arange(x-width, x+width+1),
                           np.arange(y-width, y+width+1),
                           indexing='ij')
    Xin = Xin.flatten()
    Yin = Yin.flatten()
    Xall, Yall = np.meshgrid(np.arange(np.min(Xin) - outer_width, np.max(Xin) + outer_width + 1),
                             np.arange(np.min(Yin) - outer_width, np.max(Yin) + outer_width + 1),
                             indexing='ij')
    Xall = Xall.flatten()
    Yall = Yall.flatten()

    inset = set()
    for x, y in zip(Xin, Yin):
        inset.add((x, y))

    outset = set()
    for x, y in zip(Xall, Yall):
        if (x, y) not in inset:
            outset.add((x, y))
    Xout = []
    Yout = []
    for (x, y) in outset:
        Xout.append(x)
        Yout.append(y)

    return Xout, Yout


def get_accuracy(N, si, a, b):
    """
    Standard error of the mean
    :param N: Number of photons
    :param si: standard deviation or width of the gaussian
    :param a: pixel width
    :param b: standard deviation of the background
    :return:
    """
    return np.sqrt( (si**2/N) + ((a**2/12.0)/N) + ((8*np.pi*si**4*b**2)/(a**2*N**2)))


def brightness_function(x, y, xc, yc, d, N, a, b):
    """
    Equation of the brightness
    :param x: coordinates
    :param y:
    :param xc: center of the gaussian
    :param yc:
    :param d: standard deviation
    :param N: Photon count of the roi
    :param a: pixel width
    :param b: b**2 is the bacgroun
    :return:
    """
    pi = 1.0/(2*np.pi*(d**2 + a**2/12.0))*np.exp( -1.0*a**2*((x-xc)**2 + (y-yc)**2 )/(2*(d**2 + a**2/12.0)))
    return (N*a**2)*pi + b**2


def get_gaussian_model(values):
    """
    Compute the gaussian model for a list of values
    :param values:
    :return: A dictionary
    """

    hist_values, bin_edges = np.histogram(values, bins=get_bins_number(values))
    bin_centers = [np.average(bin_edges[i:i+2]) for i, edge in enumerate(bin_edges) if i < len(bin_edges)-1]

    p0 = [np.mean(values),
          np.std(values),
          np.std(values)*np.max(hist_values)]

    popt, pcov = curve_fit(gaussian,
                           bin_centers,
                           hist_values,
                           p0=p0)

    res = {}
    res['hist_bin_centers'] = bin_centers
    res['hist_values'] = hist_values
    res['mean'], res['sdev'], res['gauss_a'] = popt
    res['interval'] = (0, np.max(values))

    return res


def get_gaussian_models(mono_model):
    """
    Compute higher gaussian models from a base model
    :param mono_model:
    :return:
    """
    N_max = int(np.power(mono_model['mean']/mono_model['sdev'], 2))

    models = [{'mean': mono_model['mean'],
               'sdev': mono_model['sdev'],
               'func': normal,
               'params': [mono_model['mean'],
                          mono_model['sdev']]}]

    for i in range(2, N_max+1):
        models.append({'mean': i * mono_model['mean'],
                       'sdev': np.sqrt(i) * mono_model['sdev'],
                       'func': normal,
                       'params': [i * mono_model['mean'],
                                  np.sqrt(i) * mono_model['sdev']]})


    return models


def normal(x, xcenter, sdev):
    return (1.0 / (sdev * np.sqrt(np.pi * 2))) *( np.exp((-1 * np.power(x - xcenter, 2)) / (2 * np.power(sdev, 2))))


def gaussian(x, xcenter, sdev, a):
    return (a / (sdev * np.sqrt(np.pi * 2))) *( np.exp((-1 * np.power(x - xcenter, 2)) / (2 * np.power(sdev, 2))))


def get_multiple_fitting(xs, ys, models, label_efficiency):
    """
    Performs a mulltiple gaussian fitting
    :param xs: x values
    :param ys: y values
    :param centers: gaussians centers
    :param sds: standard deviation of gaussians
    """

    # The area corresponding to a rectangle
    init_values = []
    for model in models:
        xc = min(xs, key=lambda x:abs(x-model['mean'])) # closest xc to the center
        yc = ys[np.where(xs == xc)][0]
        init_values.append(yc*model['sdev'])

    def mult_func(x, *coefs):
        res = 0
        for i, model in enumerate(models):
            term1 = 0
            for j in range(i, len(models)):
                term1 += binom.pmf(i+1, j+1, label_efficiency)*coefs[i]
            res += term1*model['func'](x, *model['params'])
        return res

    popt, pcov = curve_fit(mult_func, xs, ys, p0=init_values, bounds=[0, np.inf])

    res = {}
    res['params'] = popt
    res['sderrors'] = np.sqrt(np.diag(pcov))
    res['all'] = mult_func(xs, *popt)
    res['ind'] = []
    for i, model in enumerate(models):
        conv_values = popt[i]*model['func'](xs, *model['params'])
        res['ind'].append({'params': popt, 'values': conv_values})

    return res


def get_pdf_model(values):
    kde_mono = gaussian_kde(values)
    ini = 0
    end = np.max(values)
    xs_pdf_mono = np.linspace(ini, end, end+1)
    ys_pdf_mono = kde_mono.pdf(xs_pdf_mono)

    return {'pdf_function': kde_mono,
            'mean': np.mean(values),
            'sdev': np.std(values),
            'interval': (ini, end)}


def get_pdf_models(mono_model):
    N_max = int(np.power(mono_model['mean']/mono_model['sdev'], 2))
    models = []
    ini, end = mono_model['interval']
    xs = np.linspace(ini, end, end+1)
    model = {'func': mono_model['pdf_function'].pdf,
             'values': mono_model['pdf_function'].pdf(xs),
             'mean': mono_model['mean'],
             'sdev': mono_model['sdev'],
             'params': []}
    models.append(model)
    for i in range(1, N_max):
        conv_values = np.convolve(models[0]['values'], models[i-1]['values'])
        conv_xs = np.arange(len(conv_values))
        fn = interp1d(conv_xs, conv_values, fill_value=0, bounds_error=False)
        model = {'func': fn,
                 'values': conv_values,
                 'mean': (i+1)*mono_model['mean'],
                 'sdev': np.sqrt(i+1)*mono_model['sdev'],
                 'params': []}
        models.append(model)

    return models

