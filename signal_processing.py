import numpy as np
from skimage import io, draw, img_as_float
import fitting
from scipy import signal, linalg
import glob, ast


class GSignal:
    def __init__(self, parent):
        self.parent = parent
        self.hist_fixed_values = np.array([])
        self.reset()
        self.frame0_raw_brightness = 0
        self.frame0_local_bckg = 0
        self.frame0_local_bckg_avg = 0

    def generate_frame0_values(self):
        self.frame0_raw_brightness = self.get_raw_brightness(0)
        self.frame0_local_bckg_avg, self.frame0_local_bckg = self.get_local_bckg(0)

    def get_calibrated_brightness(self, bckg_method):
        """
        Brightness - background
        :param bckg_method:
        :return:
        """
        return 0

    def get_raw_brightness(self, frame):
        return 0

    def get_local_bckg(self, frame):
        return 0, 0

    def reset(self):
        self.base_values = np.array([])
        self.raw_values = np.array([])
        self.filtered_values = np.array([])
        self.kde_values = np.array([])
        self.kde_xs = np.array([])
        self.kde_ys = np.array([])
        self.kde_peak_values = np.array([])
        self.kde_steps = -1
        self.kde_fixed_values = np.array([])
        self.kde_fixed_steps = -1
        self.hist_values = np.array([])
        self.hist_steps = -1
        self.hist_fixed_values = np.array([])
        self.hist_fixed_steps = -1
        self.manual_values = np.array([])
        self.manual_steps = -1
        self.final_values = np.array([])
        self.steps_count = -1
        self.step_heights = []
        self.step_points = []

    def generate_signals(self):
        self.generate_sequence()
        if self.parent.gui_values['filter']:
            self.filter(self.parent.gui_values['median_offset'])
        if self.parent.gui_values['hist']:
            self.hist(self.parent.gui_values['bins'])
        if self.parent.gui_values['hist_fix']:
            self.hist_fix(self.parent.gui_values['bins'])
        if self.parent.gui_values['kde']:
            self.kde()
        if self.parent.gui_values['kde_fix']:
            self.kde_fix()
        if self.parent.gui_values['manual_split']:
            self.manual(self.parent.gui_values['split_points'])

        # Deciding which is final
        stepped_signals = [self.hist_values,
                           self.kde_values,
                           self.hist_fixed_values,
                           self.kde_fixed_values,
                           self.manual_values]
        for s in stepped_signals:
            if s.any():
                self.final_values = s

        self.steps_count, self.step_heights = SteppedSignal.classify(self.final_values)
        self.step_points = SteppedSignal.get_step_points(self.final_values)

    def generate_sequence(self):
        self.reset()

    def filter(self, offset):
        self.filtered_values = np.array([np.median(self.raw_values[max(0, k - offset):min(len(self.raw_values), k + offset + 1)]) for k in range(len(self.raw_values))])
        self.base_values = self.filtered_values

    def kde(self):
        """
        Compute the kde of 'values', and choose the peaks as candidates steps
        :param values:
        :return:
        """
        self.kde_xs, self.kde_ys = fitting.get_kde_values(self.base_values)
        kde_xs_peaks = signal.argrelmax(self.kde_ys, mode='wrap')[0]
        self.kde_peak_values = [self.kde_xs[i] for i in kde_xs_peaks]
        self.kde_values = SteppedSignal.get_by_step_values(self.base_values, self.kde_peak_values)
        self.kde_steps, self.kde_heights = SteppedSignal.classify(self.kde_values)

    def kde_fix(self):
        remove = False # Reset kde_value if is not active
        if not self.kde_values.any():
            self.kde()
            remove = True
        self.kde_fixed_values = SteppedSignal.get_fixed_function(self.kde_values)
        self.kde_fixed_steps, self.kde_fixed_heights = SteppedSignal.classify(self.kde_fixed_values)
        if remove:
            self.kde_values = np.array([])

    def hist(self, bins):
        """
        Generate the hist signal
        :param values:
        :param bins:
        :return:
        """
        self.hist_bin_values, self.hist_bin_edges = np.histogram(self.base_values, bins=bins, normed=True)
        self.hist_peaks = signal.argrelmax(self.hist_bin_values, mode='wrap')[0]
        if len(self.hist_peaks) == 0:
            self.hist_peaks = [np.argmax(self.hist_bin_values)]
        self.hist_peaks_values = [np.mean([self.hist_bin_edges[p], self.hist_bin_edges[p+1]]) for p in self.hist_peaks]
        self.hist_values = SteppedSignal.get_by_step_values(self.base_values, self.hist_peaks_values)
        self.hist_steps, self.hist_heights = SteppedSignal.classify(self.hist_values)

    def hist_fix(self, bins):
        remove = False # Reset hist_value if is not active
        if not self.hist_values.any():
            self.hist(bins)
            remove = True
        self.hist_fixed_values = SteppedSignal.get_fixed_function(self.hist_values)
        self.hist_fixed_steps, self.hist_fixed_heights = SteppedSignal.classify(self.hist_fixed_values)
        if remove:
            self.hist_values = np.array([])

    def manual(self, split_points):
        splits_str = '[' + str(split_points) + ']'
        split_points = ast.literal_eval(splits_str)
        self.manual_values = SteppedSignal.get_by_step_points(self.raw_values, split_points)
        self.manual_steps, self.manual_heights = SteppedSignal.classify(self.manual_values)

    def get_frame0_calibrated_brightness(self, bckg_method):
        pass


class SumSignal(GSignal):

    def generate_sequence(self):
        """
        Return for each frame of 'video' all pixel intensities within the particle region
        :param video:
        :return:
        """
        self.reset()
        intensities = []
        rr, cc = draw.disk((self.parent.x, self.parent.y), self.parent.r)
        for frame in self.parent.parent.video:
            intensities.append(sum(frame[rr, cc]))
        self.raw_values = np.array(intensities)
        self.base_values = self.raw_values

    def get_raw_brightness(self, frame):
        rr, cc = draw.disk((self.parent.x, self.parent.y), self.parent.r)
        #print("video ",self.parent.parent.video[frame][rr, cc])
        return sum(self.parent.parent.video[frame][rr, cc])

    def get_local_bckg(self, frame):
        """
        Return the intensity (by kde) of a ring of width 'width' around the particle, in frame 'frame'
        :param width:
        :return:
        """
        rr1, cc1 = draw.disk((self.parent.x, self.parent.y), self.parent.r) # inner disk
        rr2, cc2 = draw.disk((self.parent.x, self.parent.y), self.parent.r+2) # outter disk
        #Todo: set the width of the outer ring as a parameter
        s1 = set()
        for i in range(len(rr1)):
            s1.add((rr1[i], cc1[i]))
        s2 = set()
        for i in range(len(rr2)):
            s2.add((rr2[i], cc2[i]))
        intensities = []
        ring_points = s2 - s1 # only background ring points
        for p in ring_points:
            try:
                intensities.append(self.parent.parent.video[frame][p[0],p[1]])
            except IndexError: # Border was on the edge of the image
                pass
        try:
            bckg_avg = fitting.get_peak_by_kde(np.array(intensities))[0][0]
            bckg = bckg_avg*len(rr1)
            #print("background: ", bckg_avg, bckg)
            return bckg_avg, bckg
        except np.linalg.LinAlgError:
            print("background could not be estimated")
            return 0,0

    def get_frame0_calibrated_brightness(self, bckg_method):
        """
        Brightness - background
        :param bckg_method:
        :return:
        """
        if bckg_method == 'local':
            bckg = self.frame0_local_bckg
        elif bckg_method == 'global':
            rr, cc = draw.disk((self.parent.x, self.parent.y), self.parent.r)
            bckg = self.parent.parent.summary['global_bckg_avg']*(len(rr))
        return self.frame0_raw_brightness - bckg


class SteppedSignal:
    @staticmethod
    def get_by_step_values(seq, values):
        """
        Return a stepped function, for each value in 'values' return the closest value in peaks
        :param seq: sequence of all values
        :param values: reference values
        :return:
        """
        res = []
        for v in seq:
            res.append(min(values, key=lambda x: abs(x-v)))
        # Using medfilt to clean possible outliers
        return signal.medfilt(np.array(res).reshape((len(res))), 3)

    @staticmethod
    def get_fixed_function(seq):
        """
        Trying to fix up steps, keeping the longer step
        Ex:
        --__-_ -> --___
        :param seq:
        :return:
        """
        fixed_values = 0
        fix = False
        fixed_seq = np.copy(seq)
        curr = [seq[0], 0, 0] # value, range(ini, end)
        prev = None # value, range(ini, end)
        for i in range(1, len(seq)):
            if seq[i] != seq[i-1]:
                curr[2] = i-1
                if fix and prev is not None:
                    if curr[2]-curr[1] > prev[2] - prev[1]: # which step is longer
                        fixed_seq[int(prev[1]):int(prev[2]+1)] = curr[0]
                        fixed_values += (prev[2]+1) - prev[1]
                    elif curr[2]-curr[1] < prev[2] - prev[1]:
                        fixed_seq[int(curr[1]):int(curr[2]+1)] = prev[0]
                        fixed_values += (curr[2]+1) - curr[1]
                    else:
                        fixed_seq[int(curr[1]):int(curr[2]+1)] = prev[0] # if are equal always fix the upper
                        fixed_values += (curr[2]+1) - curr[1]
                prev = np.copy(curr)
                curr[0] = seq[i]
                curr[1] = i
                if seq[i] > seq[i-1]:
                    fix = True
                else:
                    fix = False
        return fixed_seq

    @staticmethod
    def get_by_step_points(seq, split_points):
        """
        Compute a stepped function from 'seq' and 'splits'.
        :param seq: sequence
        :param splits: split points of the 'value'sequence
        :return:
        """
        stepped = np.zeros(len(seq))
        new_splits = np.copy(split_points)
        new_splits = np.append(new_splits, [0, len(seq)])
        new_splits.sort()
        for i in range(1, len(new_splits)):
            stepped[new_splits[i-1]:new_splits[i]] = fitting.get_peak_by_kde(seq[new_splits[i-1]:new_splits[i]])[0][0]
        return stepped

    @staticmethod
    def get_step_points(seq):
        """
        Return the split points from a stepped function
        :param values:
        :return:
        Split points
        """
        splits = []
        for i in range(1, len(seq)):
            if seq[i] != seq[i-1]:
                splits.append(i)
        return splits

    @staticmethod
    def classify(seq):
        """
        Classify a sequence according how many downsteps it has.
        :param seq:
        :return:
        -1: monotony not decreasing (there are at least one upstep)
        0: no steps
        1: 1 step
        ...
        """
        ups = 0
        downs = 0
        steps_height = []
        for i in range(1, len(seq)):
            if seq[i] > seq[i-1]:
                ups += 1
            if seq[i] < seq[i-1]:
                downs += 1
                steps_height.append(seq[i - 1] - seq[i])
        if ups > 0:
            return -1, []
        return downs, steps_height

