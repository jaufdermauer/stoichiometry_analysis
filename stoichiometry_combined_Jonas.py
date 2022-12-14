"""
Classes and utilities for Stoichiometry Analysis
"""

__author__ = "Yuri Quintana"
__copyright__ = "Copyright 2016"
__credits__ = ["Yuri Quintana"]
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Yuri Quintana"
__email__ = "yuriqp@gmail.com"
__status__ = "Testing"

from skimage import io, draw, img_as_float
from skimage.color import rgb2gray
from skimage.morphology import reconstruction
from skimage.feature import blob_dog
from os.path import basename
from scipy import signal, linalg
from scipy.stats import gaussian_kde, binom
from scipy.optimize import curve_fit, minimize
from scipy.special import erfinv
from scipy.interpolate import interp1d
import numpy as np
import glob, ast
import math
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib import cm
from os.path import basename, relpath, realpath

import signal_processing
import fitting
import GUI
import wx

#from stoichiometry_core import StchAnalysis, StchSequence, Particle, GSignal, get_bins_number, gaussian

if __name__ == "__main__":
    app = wx.App()
    GUI.MainFrame().Show()
    app.MainLoop()



