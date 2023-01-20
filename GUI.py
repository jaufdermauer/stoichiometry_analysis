import json
import os
import wx
import wx.grid
import wx.lib.scrolledpanel as scrolled
from os.path import basename, relpath, realpath
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
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
import glob
import ast
import math
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
matplotlib.use('WXAgg')
import re
from natsort import natsorted
import csv

import pandas as pd

import analysis
import fitting

class GWidget:
    def __init__(self, widget, key, default):
        self.widget = widget
        self.default = default
        self.key = key

    def setValue(self, v):
        if isinstance(self.widget, wx.SpinCtrl):
            v = int(v)
        elif isinstance(self.widget, wx.CheckBox):
            v = v
        else:
            v = str(v)
        self.widget.SetValue(v)

    def getValue(self):
        # Todo: handle continuos updating of textctrl fields, that raises Value exceptions
        try:
            return type(self.default)(self.widget.GetValue())
        except:
            pass


class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(
            self, None, title=" Stoichiometry Analysis (Beta Version)")
        self.stch_analysis = None
        self.current_sequence = None
        self.current_particle = None
        self.current_tab = 0

        self.widgets = {}

        self.methods_key = {'sum': 'Pxls sum',
                            '2dgintegral': '2D Gauss. Integral',
                            '2dgheight': '2D Gauss. Height',
                            'local': 'Local',
                            'global': 'Global',
                            'gaussfit': 'Multiple Gaussians',
                            'pdffit': 'Multiple PDFs'}

        self.methods_val = {v: k for k, v in self.methods_key.items()}

        self.generate_menu_bar()

        p = wx.Panel(self)
        self.nb = wx.Notebook(p)

        self.generate_settings_page()
        #self.lif_to_tif()
        self.generate_detection_page()
        self.generate_particle_page()
        self.generate_photobleaching_page()
        self.generate_brightness_page()
        self.generate_stats_page()
        self.generate_analysis_page()
        self.status_bar = self.CreateStatusBar()

        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        p.SetSizer(sizer)
        self.Maximize(True)
        self.p = p
        self.stoich_calibration = -1

        #data storage in dictionary
        self.fieldnames = ["experiment","time","position","index","intensity"]
        #dictionary for characterization of spots
        self.data = {
            self.fieldnames[0] : [],    #experiment
            self.fieldnames[1] : [],    #position
            self.fieldnames[2] : [],    #time
            self.fieldnames[3] : [],     #index
            self.fieldnames[4] : []     #intensity
            }
        self.analyzed_data = {
            "time" : [],    #time
            "stoichiometry" : [],    #stoichiometry
            "stoichiometry_err" : [],
            "nparticles" : [],  #spots
            "nparticles_err" : []    
            }
        self.calibrations = []
# Pages

    def generate_menu_bar(self):

        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fitem = fileMenu.Append(wx.ID_OPEN, 'Open project')
        sitem = fileMenu.Append(wx.ID_SAVE, 'Save project')
        citem = fileMenu.Append(wx.ID_CLOSE, 'Close project')
        menubar.Append(fileMenu, '&File')
        helpMenu = wx.Menu()
        aitem = helpMenu.Append(wx.ID_ANY, 'About')
        menubar.Append(helpMenu, '&Help')
        self.SetMenuBar(menubar)
        self.Bind(wx.EVT_MENU, self.onSave, sitem)
        self.Bind(wx.EVT_MENU, self.onLoad, fitem)
        self.Bind(wx.EVT_MENU, self.onClose, citem)
        self.Bind(wx.EVT_MENU, self.onAbout, aitem)

    def onCheckbox_lts(*var):
            var = not var
            print(var)

    def generate_settings_page(self):
        self.page_settings = wx.Panel(self.nb)
        page_settings_main_box = wx.BoxSizer(wx.HORIZONTAL)

         #checkbox to save tif files in separate folder for better oversight
        cb_lts = wx.CheckBox(self.page_settings, label="save tif files to separate folders")
        cb_lts.Bind(wx.EVT_CHECKBOX, self.onCheckbox_lts(), id=cb_lts.GetId())
        cb_lts.SetValue(True)

        left_box = wx.BoxSizer(wx.VERTICAL)
        self.folder_box_ltf = wx.StaticBox(
            self.page_settings, 0, " Raw lif data to tif ")
        folder_sizer_ltf = wx.StaticBoxSizer(self.folder_box_ltf, wx.VERTICAL)
        self.folder_path_ltf = wx.TextCtrl(self.page_settings, size=(400, -1))
        folder_sizer_ltf.Add(self.folder_path_ltf, 0)
        dir_btn_ltf = wx.Button(self.page_settings, label='Browse')
        dir_btn_ltf.Bind(wx.EVT_BUTTON, lambda event: self.onSelectDataFolder_ltf(event, cb_lts.GetValue()),
                     id=dir_btn_ltf.GetId())
        folder_sizer_ltf.Add(dir_btn_ltf, 0, wx.ALIGN_RIGHT)
        
        folder_sizer_ltf.Add(cb_lts, 0, wx.ALIGN_LEFT)
        left_box.Add(folder_sizer_ltf, 0)
        left_box.AddSpacer(10)
        
        self.folder_box_cal = wx.StaticBox(
            self.page_settings, 0, " Calibration folder ")
        folder_sizer_cal = wx.StaticBoxSizer(self.folder_box_cal, wx.VERTICAL)
        self.folder_path_cal = wx.TextCtrl(self.page_settings, size=(400, -1))
        folder_sizer_cal.Add(self.folder_path_cal, 0)
        dir_btn = wx.Button(self.page_settings, label='Browse')
        dir_btn.Bind(wx.EVT_BUTTON, self.onSelectCalibrationFolder,
                     id=dir_btn.GetId())
        folder_sizer_cal.Add(dir_btn, 0, wx.ALIGN_RIGHT)

        left_box.Add(folder_sizer_cal, 0)
        left_box.AddSpacer(10)

        cb_subfolder = wx.CheckBox(self.page_settings, label="Analyze subfolders")
        cb_subfolder.Bind(wx.EVT_CHECKBOX, self.onCheckbox_lts(), id=cb_subfolder.GetId())
        cb_subfolder.SetValue(True)

        self.folder_box = wx.StaticBox(
            self.page_settings, 0, " Project's data folder ")
        folder_sizer = wx.StaticBoxSizer(self.folder_box, wx.VERTICAL)
        self.folder_path = wx.TextCtrl(self.page_settings, size=(400, -1))
        folder_sizer.Add(self.folder_path, 0)
        dir_btn = wx.Button(self.page_settings, label='Browse')
        dir_btn.Bind(wx.EVT_BUTTON, lambda event: self.onSelectDataFolder(event, cb_subfolder.GetValue()), id=dir_btn.GetId())
        folder_sizer.Add(dir_btn, 0, wx.ALIGN_RIGHT)

        folder_sizer.Add(cb_subfolder, 0, wx.ALIGN_LEFT)

        left_box.Add(folder_sizer, 0)
        left_box.AddSpacer(10)

        start_btn = wx.Button(self.page_settings, label='Start analysis!')
        start_btn.Bind(wx.EVT_BUTTON, self.onStart,
                     id=start_btn.GetId())

        general_settings_box = wx.StaticBox(
            self.page_settings, 0, " General Settings ")
        general_settings_sizer = wx.StaticBoxSizer(
            general_settings_box, wx.VERTICAL)
        gs_grid_sizer = wx.FlexGridSizer(4, 2, 2, 2)
        gs_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" Photon count coef. "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.photon_coef = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.photon_coef.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.photon_coef.GetId())
        gs_grid_sizer.Add(self.photon_coef, 0)
        gs_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" Frame rate (milliseconds)"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.frame_rate = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.frame_rate.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.frame_rate.GetId())
        gs_grid_sizer.Add(self.frame_rate, 0)
        gs_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" Pixel size (nanometers) "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.pixel_size = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.pixel_size.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.pixel_size.GetId())
        gs_grid_sizer.Add(self.pixel_size, 0)
        gs_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" Particle ROI radius (pixels)"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.roi_radius = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.roi_radius.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.roi_radius.GetId())
        gs_grid_sizer.Add(self.roi_radius, 0)
        general_settings_sizer.Add(gs_grid_sizer, 0)
        left_box.Add(general_settings_sizer, 0, wx.EXPAND)
        left_box.AddSpacer(10)

        defaults_det_grid_sizer = wx.FlexGridSizer(4, 2, 2, 2)
        defaults_det_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" min_sigma "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.default_min_sigma = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.default_min_sigma.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.default_min_sigma.GetId())
        defaults_det_grid_sizer.Add(self.default_min_sigma, 0)
        defaults_det_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" max_sigma "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.default_max_sigma = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.default_max_sigma.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.default_max_sigma.GetId())
        defaults_det_grid_sizer.Add(self.default_max_sigma, 0)

        defaults_det_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" threshold "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.default_threshold = wx.TextCtrl(self.page_settings, size=(80, -1))
        self.default_threshold.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.default_threshold.GetId())
        defaults_det_grid_sizer.Add(self.default_threshold, 0, wx.EXPAND)

        defaults_det_grid_sizer.Add(wx.StaticText(
            self.page_settings, label=" iterations "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.iterations = wx.TextCtrl(
            self.page_settings, size=(80, -1), value='4')
        self.iterations.Bind(
            wx.EVT_TEXT, self.onChangeSettingsField, id=self.iterations.GetId())
        defaults_det_grid_sizer.Add(self.iterations, 0)

        part_detection_box = wx.StaticBox(
            self.page_settings, 0, " Particles Detection ")
        part_detection_sizer = wx.StaticBoxSizer(
            part_detection_box, wx.VERTICAL)
        part_detection_sizer.Add(defaults_det_grid_sizer, 0, wx.EXPAND)

        defaults_pro_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        defaults_pro_grid_sizer.Add(wx.StaticText(
            self.page_settings, label="Median Offset"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.default_offset = wx.ComboBox(self.page_settings, choices=[
                                          '1', '3', '5', '7', '9', '11'], size=(80, -1))
        self.default_offset.Bind(
            wx.EVT_COMBOBOX, self.onChangeSettingsField, id=self.default_offset.GetId())
        defaults_pro_grid_sizer.Add(self.default_offset, 0)
        defaults_pro_grid_sizer.Add(wx.StaticText(
            self.page_settings, label="Bins"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.default_bins = wx.SpinCtrl(
            self.page_settings, min=1, max=120, size=(80, -1))
        self.default_bins.Bind(
            wx.EVT_SPINCTRL, self.onChangeSettingsField, id=self.default_bins.GetId())
        defaults_pro_grid_sizer.Add(self.default_bins, 0)

        part_processing_box = wx.StaticBox(
            self.page_settings, 0, " Particle Processing ")
        part_processing_sizer = wx.StaticBoxSizer(
            part_processing_box, wx.VERTICAL)
        part_processing_sizer.Add(defaults_pro_grid_sizer)

        defaults_box = wx.StaticBox(
            self.page_settings, 0, " Project defaults ")
        defaults_sizer = wx.StaticBoxSizer(defaults_box, wx.VERTICAL)
        defaults_sizer.Add(part_detection_sizer, 0, wx.EXPAND)
        defaults_sizer.Add(part_processing_sizer, 0, wx.EXPAND)
        left_box.Add(defaults_sizer, 0, wx.EXPAND)
        page_settings_main_box.Add(left_box, 0)

        text = " Project's files (check files to be used for calibration)"
        self.setting_box = wx.StaticBox(self.page_settings, 1, text)
        self.settings_box_sizer = wx.StaticBoxSizer(
            self.setting_box, wx.VERTICAL)
        self.scrolled_panel = scrolled.ScrolledPanel(self.page_settings)
        self.scrolled_panel.SetAutoLayout(1)
        self.scrolled_panel.SetupScrolling()
        self.scrolled_panel_sizer = wx.BoxSizer(wx.VERTICAL)
        self.scrolled_panel.SetSizer(self.scrolled_panel_sizer)
        self.settings_box_sizer.Add(self.scrolled_panel, 1, wx.EXPAND)
        page_settings_main_box.Add(self.settings_box_sizer, 1, wx.EXPAND)

        self.page_settings.SetSizer(page_settings_main_box, wx.EXPAND)
        self.nb.AddPage(self.page_settings, "Project Settings")

        defaults = [('photon_coef', self.photon_coef, 12.5),
                    ('frame_rate', self.frame_rate, 60),
                    ('pixel_size', self.pixel_size, 100),
                    ('def_roi_radius', self.roi_radius, 5),
                    ('def_min_sigma', self.default_min_sigma, 2.0),
                    ('def_max_sigma', self.default_max_sigma, 3.0),
                    ('def_threshold', self.default_threshold, 0.000025),
                    ('def_iterations', self.iterations, 4),
                    ('def_median_offset', self.default_offset, 3),
                    ('def_bins', self.default_bins, 10)]

        self.analysis_widgets = []

        for key, widget, value in defaults:
            gwidget = GWidget(widget, key, value)
            gwidget.setValue(value)
            self.widgets[widget.GetId()] = gwidget
            self.analysis_widgets.append(gwidget)

        # start button
        left_box.AddSpacer(10)
        left_box.Add(start_btn, 0)
        left_box.AddSpacer(10)

    def generate_detection_page(self):

        self.page_detection = wx.Panel(self.nb)
        self.nb.AddPage(self.page_detection, "Particle Detection")
        page_detection_main_box = wx.BoxSizer(wx.HORIZONTAL)

        self.page_detection_image_sizer = wx.BoxSizer(wx.VERTICAL)
        self.page_detection_image_panel = wx.Panel(self.page_detection)
        self.page_detection_image_sizer.Add(
            self.page_detection_image_panel, 1, wx.EXPAND)
        page_detection_main_box.Add(
            self.page_detection_image_sizer, 1, wx.EXPAND)

        menu_sizer = wx.BoxSizer(wx.VERTICAL)
        menu_sizer.AddSpacer(5)
        # Todo: scroll panel

        replicate_box = wx.StaticBox(self.page_detection, 0, " Replicate ")
        replicate_sizer = wx.StaticBoxSizer(replicate_box, wx.VERTICAL)
        replicate_sizer.AddSpacer(5)

        self.replicate_combo = wx.ComboBox(
            self.page_detection, choices=[], size=(350, -1))
        self.replicate_combo.Bind(
            wx.EVT_COMBOBOX, self.onChangeImage, id=self.replicate_combo.GetId())
        replicate_sizer.Add(self.replicate_combo, 0)
        replicate_sizer.AddSpacer(10)

        frame_sizer = wx.BoxSizer(wx.HORIZONTAL)
        frame_sizer.Add(wx.StaticText(self.page_detection,
                        label=" Frame "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.frame = wx.SpinCtrl(
            self.page_detection, value='0', min=0, max=0, size=(80, -1))
        frame_sizer.Add(self.frame, 1)
        self.frame.Bind(wx.EVT_SPINCTRL, self.onChangeFrame,
                        id=self.frame.GetId())
        frame_sizer.AddSpacer(15)
        frame_sizer.Add(wx.StaticText(self.page_detection,
                        label=' Calibration: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.calibration = wx.StaticText(self.page_detection, label='-')
        frame_sizer.Add(self.calibration, 1, wx.ALIGN_CENTER_VERTICAL)
        replicate_sizer.Add(frame_sizer, 1)
        replicate_sizer.AddSpacer(10)

        detection_box = wx.StaticBox(
            self.page_detection, 0, " Detection and localization ")
        detection_sizer = wx.StaticBoxSizer(detection_box, wx.VERTICAL)
        detection_sizer.AddSpacer(5)

        detection_line_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.dog_box = wx.StaticBox(self.page_detection, 0, " Automatic ")
        dog_sizer = wx.StaticBoxSizer(self.dog_box, wx.VERTICAL)

        dog_params_sizer = wx.FlexGridSizer(3, 2, 2, 2)
        dog_params_sizer.Add(wx.StaticText(
            self.page_detection, label=" min_sigma "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.min_sigma = wx.TextCtrl(
            self.page_detection, value='2.0', size=(80, -1))
        self.min_sigma.Bind(
            wx.EVT_TEXT, self.onChangeDetectionField, id=self.min_sigma.GetId())
        dog_params_sizer.Add(self.min_sigma, 0)
        dog_params_sizer.Add(wx.StaticText(
            self.page_detection, label=" max_sigma "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.max_sigma = wx.TextCtrl(
            self.page_detection, value='3.0', size=(80, -1))
        self.max_sigma.Bind(
            wx.EVT_TEXT, self.onChangeDetectionField, id=self.max_sigma.GetId())
        dog_params_sizer.Add(self.max_sigma, 0)

        dog_params_sizer.Add(wx.StaticText(
            self.page_detection, label=" threshold "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.threshold = wx.TextCtrl(
            self.page_detection, size=(80, -1), value='0.0001')
        self.threshold.Bind(
            wx.EVT_TEXT, self.onChangeDetectionField, id=self.threshold.GetId())
        dog_params_sizer.Add(self.threshold, 0)
        
        dog_sizer.Add(dog_params_sizer, 0)
        
        dog_sizer.AddSpacer(5)
        self.detect_btn = wx.Button(self.page_detection, label='Detect (DoG)')
        self.detect_btn.Bind(wx.EVT_BUTTON, self.onDetect,
                             id=self.detect_btn.GetId())
        dog_sizer.Add(self.detect_btn, 0, wx.ALIGN_RIGHT)
        dog_sizer.Add(wx.StaticLine(self.page_detection),
                      0, wx.ALL | wx.EXPAND, 5)
        self.background_btn = wx.Button(
            self.page_detection, label='Background')
        self.background_btn.Bind(
            wx.EVT_BUTTON, self.onBackground, id=self.background_btn.GetId())
        dog_sizer.Add(self.background_btn, 0, wx.ALIGN_RIGHT)

        detection_line_sizer.Add(dog_sizer, 1)

        manual_roi_box = wx.StaticBox(self.page_detection, 0, " ROI handling ")
        manual_roi_sizer = wx.StaticBoxSizer(manual_roi_box, wx.VERTICAL)

        show_rois_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.show_check_box = wx.CheckBox(
            self.page_detection, label='Show ROIs ')
        self.show_check_box.Bind(
            wx.EVT_CHECKBOX, self.onShowParticles, id=self.show_check_box.GetId())
        show_rois_sizer.Add(self.show_check_box, 0)
        manual_roi_sizer.Add(show_rois_sizer, 0)
        manual_roi_sizer.Add(wx.StaticLine(
            self.page_detection), 0, wx.ALL | wx.EXPAND, 5)

        self.manual_roi_activation = wx.CheckBox(
            self.page_detection, label='Activate')
        self.manual_roi_activation.Bind(wx.EVT_CHECKBOX, self.onManualParticleHandling,
                                        id=self.manual_roi_activation.GetId())
        manual_roi_sizer.Add(self.manual_roi_activation, 0)
        self.remove_opt = wx.RadioButton(
            self.page_detection, label='Remove', style=wx.RB_GROUP)
        manual_roi_sizer.Add(self.remove_opt, 0)
        self.add_opt = wx.RadioButton(self.page_detection, label='Add')
        manual_roi_sizer.Add(self.add_opt, 0)
        self.add_type = wx.ComboBox(self.page_detection, value='Particle', choices=[
                                    'Particle', 'Background'], size=((100, -1)))
        self.add_radius = wx.SpinCtrl(
            self.page_detection, value='3', min=1, max=40, size=(50, -1))
        add_opt_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        add_opt_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Type"), 0)
        add_opt_grid_sizer.Add(self.add_type, 0)
        add_opt_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Width"), 0)
        add_opt_grid_sizer.Add(self.add_radius, 0)
        manual_roi_sizer.Add(add_opt_grid_sizer, 0, wx.ALIGN_RIGHT)
        detection_line_sizer.Add(manual_roi_sizer, 1)

        self.manual_roi_controls = [self.remove_opt,
                                    self.add_opt,
                                    self.add_radius,
                                    self.add_type]

        detection_sizer.Add(detection_line_sizer, 0, wx.EXPAND)
        detection_sizer.AddSpacer(10)

        discard_box = wx.StaticBox(self.page_detection, 0, " Discard ")
        discard_sizer = wx.StaticBoxSizer(discard_box, wx.VERTICAL)

        dist_discard = wx.BoxSizer(wx.HORIZONTAL)
        dist_discard.Add(wx.StaticText(self.page_detection,
                         label='Distance: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.discard_btn = wx.Button(self.page_detection, label='Discard')
        self.discard_btn.Bind(wx.EVT_BUTTON, self.onDiscard,
                              id=self.discard_btn.GetId())
        dist_discard.Add(self.discard_btn, 1)

        discard_sizer.Add(dist_discard, 0)
        discard_sizer.Add(wx.StaticLine(self.page_detection),
                          0, wx.ALL | wx.EXPAND, 5)

        discard_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        discard_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Threshold: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.width_threshold = wx.TextCtrl(
            self.page_detection, size=((80, -1)), value='0.95')
        discard_grid_sizer.Add(self.width_threshold, 1)
        self.discard_preview = wx.CheckBox(
            self.page_detection, label='Preview')
        self.discard_preview.Bind(
            wx.EVT_CHECKBOX, self.onDiscardPreview, id=self.discard_preview.GetId())
        discard_grid_sizer.Add(self.discard_preview, 0,
                               wx.ALIGN_CENTER_VERTICAL)
        self.wdiscard_btn = wx.Button(
            self.page_detection, label='Discard', size=((80, -1)))
        self.wdiscard_btn.Bind(
            wx.EVT_BUTTON, self.onWDiscard, id=self.wdiscard_btn.GetId())
        discard_grid_sizer.Add(self.wdiscard_btn, 1)

        discard_sizer.Add(discard_grid_sizer, 0)

        discard_sizer.Add(wx.StaticLine(self.page_detection),
                          0, wx.ALL | wx.EXPAND, 5)

        stack_discard = wx.BoxSizer(wx.HORIZONTAL)
        stack_discard.Add(wx.StaticText(self.page_detection,
                          label='Stack: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sdiscard_btn = wx.Button(self.page_detection, label='Discard')
        self.sdiscard_btn.Bind(
            wx.EVT_BUTTON, self.onStackDiscard, id=self.sdiscard_btn.GetId())
        stack_discard.Add(self.sdiscard_btn, 1)

        discard_sizer.Add(stack_discard, 0)

        localize_box = wx.StaticBox(self.page_detection, 0, " Localization ")
        localize_sizer = wx.StaticBoxSizer(localize_box, wx.HORIZONTAL)

        localize_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        localize_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Method: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.fit_method = wx.ComboBox(self.page_detection, choices=[
                                      'LS', 'MLE'], value='LS')
        localize_grid_sizer.Add(self.fit_method, 0)
        localize_grid_sizer.Add(wx.StaticText(self.page_detection), 0)
        self.localize_btn = wx.Button(self.page_detection, label='Localize')
        self.localize_btn.Bind(
            wx.EVT_BUTTON, self.onLocalize, id=self.localize_btn.GetId())
        localize_grid_sizer.Add(self.localize_btn, 0)
        localize_sizer.Add(localize_grid_sizer, 1)

        discard_localize_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        detection_sizer.Add(discard_localize_line_sizer, 0, wx.EXPAND)
        discard_localize_line_sizer.Add(discard_sizer, 1)
        discard_localize_line_sizer.Add(localize_sizer, 1)

        Iterate_box = wx.StaticBox(self.page_detection, 0, " Iteration ")
        Iterate_sizer = wx.StaticBoxSizer(Iterate_box, wx.HORIZONTAL)

        Iterate_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        Iterate_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Method: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.fit_method = wx.ComboBox(self.page_detection, choices=[
                                      'LS', 'MLE'], value='LS')
        Iterate_grid_sizer.Add(self.fit_method, 0)
        Iterate_grid_sizer.Add(wx.StaticText(self.page_detection), 0)
        self.Iterate_btn = wx.Button(self.page_detection, label='Iterate')
        self.Iterate_btn.Bind(wx.EVT_BUTTON, self.onIterate,
                              id=self.Iterate_btn.GetId())
        Iterate_grid_sizer.Add(self.Iterate_btn, 0)
        Iterate_sizer.Add(Iterate_grid_sizer, 1)

        discard_Iterate_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        detection_sizer.Add(discard_Iterate_line_sizer, 0, wx.EXPAND)
        discard_Iterate_line_sizer.Add(discard_sizer, 1)
        discard_Iterate_line_sizer.Add(Iterate_sizer, 1)

        self.summary_box = wx.StaticBox(self.page_detection, 0, " Summary (0)")
        summary_box_sizer = wx.StaticBoxSizer(self.summary_box, wx.VERTICAL)

        fitted_width_box = wx.StaticBox(
            self.page_detection, 0, " Fitted width ")
        fitted_width_sizer = wx.StaticBoxSizer(fitted_width_box, wx.HORIZONTAL)
        fitted_width_grid_sizer = wx.FlexGridSizer(1, 4, 2, 2)
        fitted_width_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Mean: "), 1, wx.EXPAND)
        self.fitted_width_avg = wx.StaticText(self.page_detection, label='-')
        fitted_width_grid_sizer.Add(self.fitted_width_avg, 1)
        fitted_width_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Std: "), 1)
        self.fitted_width_std = wx.StaticText(self.page_detection, label='-')
        fitted_width_grid_sizer.Add(self.fitted_width_std, 1)
        fitted_width_grid_sizer.AddGrowableCol(1, 1)
        fitted_width_grid_sizer.AddGrowableCol(3, 1)
        fitted_width_sizer.Add(fitted_width_grid_sizer, 1)

        background_box = wx.StaticBox(self.page_detection, 0, " Background ")
        background_box_sizer = wx.StaticBoxSizer(background_box, wx.HORIZONTAL)

        bckg_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        bckg_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Global:"))
        self.global_bckg_avg = wx.StaticText(self.page_detection, label="-")
        bckg_grid_sizer.Add(self.global_bckg_avg, 0)
        bckg_grid_sizer.Add(wx.StaticText(self.page_detection, label="Local:"))
        self.local_bckg_avg = wx.StaticText(self.page_detection, label="-")
        bckg_grid_sizer.Add(self.local_bckg_avg, 0)
        background_box_sizer.Add(bckg_grid_sizer, 1)

        accuracy_box = wx.StaticBox(self.page_detection, 0, " Accuracy ")
        accuracy_box_sizer = wx.StaticBoxSizer(accuracy_box, wx.HORIZONTAL)

        accr_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        accr_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Global:"))
        self.accr_global_avg = wx.StaticText(self.page_detection, label="-")
        accr_grid_sizer.Add(self.accr_global_avg, 0)
        accr_grid_sizer.Add(wx.StaticText(self.page_detection, label="Local:"))
        self.accr_local_avg = wx.StaticText(self.page_detection, label="-")
        accr_grid_sizer.Add(self.accr_local_avg, 0)
        accuracy_box_sizer.Add(accr_grid_sizer, 1)

        bckg_accr_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bckg_accr_line_sizer.Add(background_box_sizer, 1)
        bckg_accr_line_sizer.Add(accuracy_box_sizer, 1)

        summary_box_sizer.AddSpacer(5)
        summary_box_sizer.Add(fitted_width_sizer, 0, wx.EXPAND)
        summary_box_sizer.AddSpacer(5)
        summary_box_sizer.Add(bckg_accr_line_sizer, 0, wx.EXPAND)

        sel_particle_box = wx.StaticBox(self.page_detection, 0, " Particle ")
        sel_particle_box_sizer = wx.StaticBoxSizer(
            sel_particle_box, wx.VERTICAL)

        part_flex_grid_sizer = wx.FlexGridSizer(1, 2, 2, 2)
        part_flex_grid_sizer.AddGrowableCol(0, 1)
        part_flex_grid_sizer.AddGrowableCol(1, 1)

        part_general_box = wx.StaticBox(self.page_detection, 0, " General ")
        part_general_box_sizer = wx.StaticBoxSizer(
            part_general_box, wx.HORIZONTAL)

        part_general_grid_sizer = wx.FlexGridSizer(4, 2, 2, 2)
        part_general_grid_sizer.Add(
            wx.StaticText(self.page_detection, label="x: "))
        self.coords_x_part = wx.StaticText(self.page_detection, label="-")
        part_general_grid_sizer.Add(self.coords_x_part, 0)
        part_general_grid_sizer.Add(
            wx.StaticText(self.page_detection, label="y: "))
        self.coords_y_part = wx.StaticText(self.page_detection, label="-")
        part_general_grid_sizer.Add(self.coords_y_part, 0)
        part_general_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Fit width: "))
        self.fitted_width = wx.StaticText(self.page_detection, label="-")
        part_general_grid_sizer.Add(self.fitted_width, 0)
        visualize_btn = wx.Button(
            self.page_detection, label='Visualize', size=((80, -1)))
        visualize_btn.Bind(wx.EVT_BUTTON, self.onVisualize,
                           id=visualize_btn.GetId())
        part_general_grid_sizer.Add(visualize_btn, 0)
        part_general_grid_sizer.Add(
            wx.StaticText(self.page_detection, label=""), 0)

        part_general_box_sizer.Add(part_general_grid_sizer, 0, wx.EXPAND)

        part_accuracy_box = wx.StaticBox(self.page_detection, 0, " Accuracy ")
        part_accuracy_box_sizer = wx.StaticBoxSizer(
            part_accuracy_box, wx.HORIZONTAL)

        part_accr_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        part_accr_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Global:"))
        self.part_accr_global = wx.StaticText(self.page_detection, label="-")
        part_accr_grid_sizer.Add(self.part_accr_global, 0)
        part_accr_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Local:"))
        self.part_accr_local = wx.StaticText(self.page_detection, label="-")
        part_accr_grid_sizer.Add(self.part_accr_local, 0)
        part_accuracy_box_sizer.Add(part_accr_grid_sizer, 0, wx.EXPAND)

        part_flex_grid_sizer.Add(part_general_box_sizer, 0, wx.EXPAND)
        part_flex_grid_sizer.Add(part_accuracy_box_sizer, 0, wx.EXPAND)
        sel_particle_box_sizer.Add(part_flex_grid_sizer, 0, wx.EXPAND)

        part_brightness_box = wx.StaticBox(
            self.page_detection, 0, " Brightness (raw) and Background (avg)")
        part_brightness_box_sizer = wx.StaticBoxSizer(
            part_brightness_box, wx.HORIZONTAL)

        part_brightness_grid_sizer = wx.FlexGridSizer(3, 3, 2, 2)
        part_brightness_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="Pixels sum"))
        self.part_brightness_sum = wx.StaticText(
            self.page_detection, label="-")
        part_brightness_grid_sizer.Add(self.part_brightness_sum, 1, wx.EXPAND)
        self.part_brightness_sum_bckg = wx.StaticText(
            self.page_detection, label="-")
        part_brightness_grid_sizer.Add(
            self.part_brightness_sum_bckg, 1, wx.EXPAND)
        part_brightness_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="2D Gauss. height"))
        self.part_brightness_2dheight = wx.StaticText(
            self.page_detection, label="-")
        part_brightness_grid_sizer.Add(
            self.part_brightness_2dheight, 1, wx.EXPAND)
        self.part_brightness_2dheight_bckg = wx.StaticText(
            self.page_detection, label="-")
        part_brightness_grid_sizer.Add(
            self.part_brightness_2dheight_bckg, 1, wx.EXPAND)
        part_brightness_grid_sizer.Add(wx.StaticText(
            self.page_detection, label="2D Gauss. integral"))
        self.part_brightness_2dintegral = wx.StaticText(
            self.page_detection, label="-")
        part_brightness_grid_sizer.Add(
            self.part_brightness_2dintegral, 1, wx.EXPAND)
        self.part_brightness_2dintegral_bckg = wx.StaticText(
            self.page_detection, label="-")
        part_brightness_grid_sizer.Add(
            self.part_brightness_2dintegral_bckg, 1, wx.EXPAND)
        part_brightness_box_sizer.Add(part_brightness_grid_sizer, 1, wx.EXPAND)

        part_brightness_grid_sizer.AddGrowableCol(1, 1)
        part_brightness_grid_sizer.AddGrowableCol(2, 1)

        sel_particle_box_sizer.Add(part_brightness_box_sizer, 1, wx.EXPAND)

        replicate_sizer.Add(detection_sizer, 0, wx.EXPAND)
        replicate_sizer.AddSpacer(5)
        replicate_sizer.Add(summary_box_sizer, 0, wx.EXPAND)
        replicate_sizer.AddSpacer(5)
        replicate_sizer.Add(sel_particle_box_sizer, 0, wx.EXPAND)
        menu_sizer.Add(replicate_sizer, 0)

        page_detection_main_box.Add(menu_sizer, 0)

        self.page_detection.SetSizer(page_detection_main_box, wx.EXPAND)

        defaults = [('frame', self.frame, 0),
                    ('min_sigma', self.min_sigma, float(
                        self.default_min_sigma.GetValue())),
                    ('max_sigma', self.max_sigma, float(
                        self.default_max_sigma.GetValue())),
                    ('threshold', self.threshold, float(
                        self.default_threshold.GetValue())),
                    ('roi_radius', self.add_radius, int(self.roi_radius.GetValue()))]

        self.sequence_widgets = []

        for key, widget, value in defaults:
            gwidget = GWidget(widget, key, value)
            gwidget.setValue(value)
            self.widgets[widget.GetId()] = gwidget
            self.sequence_widgets.append(gwidget)

    def generate_particle_page(self):
        self.page_processing = wx.Panel(self.nb)
        self.nb.AddPage(self.page_processing, "Particle Processing")
        self.nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED,
                     self.onChangeTab, id=self.nb.GetId())

        page_part_main_box = wx.BoxSizer(wx.HORIZONTAL)

        self.page_part_image_sizer = wx.BoxSizer(wx.VERTICAL)
        self.page_part_image_panel = wx.Panel(self.page_processing)
        self.page_part_image_sizer.Add(
            self.page_part_image_panel, 1, wx.EXPAND)
        page_part_main_box.Add(self.page_part_image_sizer, 1, wx.EXPAND)

        page_part_menu_sizer = wx.BoxSizer(wx.VERTICAL)
        page_part_menu_sizer.AddSpacer(10)

        particle_box = wx.StaticBox(self.page_processing, 0, "Particle:")
        particle_sizer = wx.StaticBoxSizer(particle_box, wx.HORIZONTAL)
        id_sizer = wx.FlexGridSizer(2, 2, 2, 2)

        id_sizer.Add(wx.StaticText(self.page_processing,
                     label="Id"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.particle_id = wx.SpinCtrl(
            self.page_processing, value='1', min=0, max=0, size=((80, -1)))
        self.particle_id.Bind(
            wx.EVT_SPINCTRL, self.onChangeParticle, id=self.particle_id.GetId())
        id_sizer.Add(self.particle_id, 0)
        id_sizer.Add(wx.StaticText(self.page_processing), 0)
        self.dlt_btn = wx.Button(self.page_processing,
                                 label='Delete', size=((80, -1)))
        self.dlt_btn.Bind(wx.EVT_BUTTON, self.onDeleteParticle,
                          id=self.dlt_btn.GetId())
        id_sizer.Add(self.dlt_btn, 0)

        coords_box = wx.StaticBox(self.page_processing, 0, " Info: ")
        coords_sizer = wx.StaticBoxSizer(coords_box, wx.HORIZONTAL)
        coords_grid_sizer = wx.FlexGridSizer(3, 2, 2, 4)
        coords_grid_sizer.Add(wx.StaticText(
            self.page_processing, label="x: "), 0, wx.ALIGN_RIGHT)
        self.coords_x = wx.StaticText(self.page_processing, label="-")
        coords_grid_sizer.Add(self.coords_x, 1)
        coords_grid_sizer.Add(wx.StaticText(
            self.page_processing, label="y: "), 0, wx.ALIGN_RIGHT)
        self.coords_y = wx.StaticText(self.page_processing, label="-")
        coords_grid_sizer.Add(self.coords_y, 1)
        coords_grid_sizer.Add(wx.StaticText(
            self.page_processing, label="width: "), 0, wx.ALIGN_RIGHT)
        self.width = wx.StaticText(self.page_processing, label="-")
        coords_grid_sizer.Add(self.width, 1)
        coords_sizer.Add(coords_grid_sizer)

        particle_sizer.Add(id_sizer, 0)
        particle_sizer.Add(coords_sizer, 1)

        page_part_menu_sizer.Add(particle_sizer, 0, wx.EXPAND)
        page_part_menu_sizer.AddSpacer(10)

        brightness_box = wx.StaticBox(self.page_processing, 0, " Brightness ")
        brightness_sizer = wx.StaticBoxSizer(brightness_box, wx.VERTICAL)
        bg_method_grid_sizer = wx.FlexGridSizer(1, 2, 2, 2)
        bg_method_grid_sizer.Add(wx.StaticText(
            self.page_processing, label="Method:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.bg_method = wx.ComboBox(self.page_processing, value=self.methods_key['sum'],
                                     choices=[self.methods_key['sum'],
                                              self.methods_key['2dgintegral'],
                                              self.methods_key['2dgheight']])
        bg_method_grid_sizer.Add(self.bg_method, 0)
        brightness_sizer.Add(bg_method_grid_sizer, 0, wx.EXPAND)
        # self.background = wx.CheckBox(self.page_processing, label='Background')
        # self.background.Bind(wx.EVT_CHECKBOX, self.onChangeParticleField, id=self.background.GetId())
        # brightness_sizer.Add(self.background, 0)
        page_part_menu_sizer.Add(brightness_sizer, 0, wx.EXPAND)
        page_part_menu_sizer.AddSpacer(10)

        self.filter_box = wx.StaticBox(self.page_processing, 0, " Filtering ")
        filter_sizer = wx.StaticBoxSizer(self.filter_box, wx.VERTICAL)
        offset_sizer = wx.FlexGridSizer(1, 2, 2, 2)
        offset_sizer.Add(wx.StaticText(self.page_processing,
                         label="Median Offset"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.offset = wx.ComboBox(self.page_processing, value='3', choices=[
                                  '1', '3', '5', '7', '9', '11'], size=(80, -1))
        self.offset.Bind(
            wx.EVT_COMBOBOX, self.onChangeParticleField, id=self.offset.GetId())
        offset_sizer.Add(self.offset, 0)
        filter_sizer.Add(offset_sizer, 0)
        self.median_filter = wx.CheckBox(
            self.page_processing, label="Median Filter")
        self.median_filter.Bind(
            wx.EVT_CHECKBOX, self.onChangeParticleField, id=self.median_filter.GetId())
        filter_sizer.Add(self.median_filter, 0)
        page_part_menu_sizer.Add(filter_sizer, 0, wx.EXPAND)
        page_part_menu_sizer.AddSpacer(10)

        steps_box = wx.StaticBox(self.page_processing, 0, "Step Detection:")
        steps_sizer = wx.StaticBoxSizer(steps_box, wx.VERTICAL)
        steps_sizer.AddSpacer(5)

        self.kde_box = wx.StaticBox(
            self.page_processing, 0, "Kernel Density Function")
        kde_sizer = wx.StaticBoxSizer(self.kde_box, wx.VERTICAL)
        kde_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        kde_sizer.Add(kde_grid_sizer, 0, wx.EXPAND)
        self.kde = wx.CheckBox(self.page_processing, label="KDE ")
        self.kde.Bind(wx.EVT_CHECKBOX, self.onChangeParticleField,
                      id=self.kde.GetId())
        kde_grid_sizer.Add(self.kde, 0)
        self.kde_steps_count = wx.StaticText(self.page_processing, label='-')
        kde_grid_sizer.Add(self.kde_steps_count, 0, wx.ALIGN_CENTER_VERTICAL)
        self.kde_fix = wx.CheckBox(self.page_processing, label="Fix KDE ")
        self.kde_fix.Bind(
            wx.EVT_CHECKBOX, self.onChangeParticleField, id=self.kde_fix.GetId())
        kde_grid_sizer.Add(self.kde_fix, 0)
        self.kde_fix_steps_count = wx.StaticText(
            self.page_processing, label='-')
        kde_grid_sizer.Add(self.kde_fix_steps_count,
                           0, wx.ALIGN_CENTER_VERTICAL)
        steps_sizer.Add(kde_sizer, 0, wx.EXPAND)
        steps_sizer.AddSpacer(5)

        self.hist_box = wx.StaticBox(self.page_processing, 0, "Histogram")
        hist_sizer = wx.StaticBoxSizer(self.hist_box, wx.VERTICAL)
        bins_sizer = wx.FlexGridSizer(1, 3, 2, 2)
        bins_sizer.Add(wx.StaticText(self.page_processing,
                       label="Bins"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.bins = wx.SpinCtrl(self.page_processing,
                                value='10', min=1, max=120, size=(80, -1))
        self.bins.Bind(wx.EVT_SPINCTRL,
                       self.onChangeParticleField, id=self.bins.GetId())
        bins_sizer.Add(self.bins, 0)
        self.auto_btn = wx.Button(self.page_processing, label="Auto")
        self.auto_btn.Bind(wx.EVT_BUTTON, self.onAutoBin,
                           id=self.auto_btn.GetId())
        bins_sizer.Add(self.auto_btn, 0)
        hist_sizer.Add(bins_sizer, 0)
        hist_grid_sizer = wx.FlexGridSizer(2, 2, 2, 2)
        hist_sizer.Add(hist_grid_sizer, 0, wx.EXPAND)
        self.hist = wx.CheckBox(self.page_processing, label="Histogram ")
        self.hist.Bind(wx.EVT_CHECKBOX,
                       self.onChangeParticleField, id=self.hist.GetId())
        hist_grid_sizer.Add(self.hist, 0)
        self.hist_steps_count = wx.StaticText(self.page_processing, label='-')
        hist_grid_sizer.Add(self.hist_steps_count, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hist_fix = wx.CheckBox(
            self.page_processing, label="Fix Histogram ")
        self.hist_fix.Bind(
            wx.EVT_CHECKBOX, self.onChangeParticleField, id=self.hist_fix.GetId())
        hist_grid_sizer.Add(self.hist_fix, 0)
        self.hist_fix_steps_count = wx.StaticText(
            self.page_processing, label='-')
        hist_grid_sizer.Add(self.hist_fix_steps_count,
                            0, wx.ALIGN_CENTER_VERTICAL)
        steps_sizer.Add(hist_sizer, 0, wx.EXPAND)
        steps_sizer.AddSpacer(5)

        manual_box = wx.StaticBox(self.page_processing, 0, "Manual")
        manual_sizer = wx.StaticBoxSizer(manual_box, wx.VERTICAL)
        manual_sizer.Add(wx.StaticText(self.page_processing,
                         label='Split points (, separated):'), 0)
        self.split_points = wx.TextCtrl(self.page_processing, size=(200, -1))
        self.split_points.Bind(
            wx.EVT_TEXT, self.onChangeSplitPoints, id=self.split_points.GetId())
        manual_sizer.Add(self.split_points, 0)
        manual_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        manual_sizer.Add(manual_line_sizer, 0, wx.EXPAND)
        self.manual_split = wx.CheckBox(self.page_processing, label="Manual ")
        self.manual_split.Bind(
            wx.EVT_CHECKBOX, self.onChangeParticleField, id=self.manual_split.GetId())
        manual_line_sizer.Add(self.manual_split, 0)
        self.manual_steps_count = wx.StaticText(
            self.page_processing, label='-')
        manual_line_sizer.Add(self.manual_steps_count, 0,
                              wx.ALIGN_CENTER_VERTICAL)
        steps_sizer.Add(manual_sizer, 0, wx.EXPAND)

        page_part_menu_sizer.Add(steps_sizer, 0, wx.EXPAND)

        summary_box = wx.StaticBox(self.page_processing, 0, " Summary ")
        summary_sizer = wx.StaticBoxSizer(summary_box, wx.VERTICAL)

        count_grid_sizer = wx.FlexGridSizer(1, 2, 2, 2)
        count_grid_sizer.Add(wx.StaticText(
            self.page_processing, label="Steps count: "))
        self.steps_count = wx.StaticText(self.page_processing, label="-")
        count_grid_sizer.Add(self.steps_count, 0)
        summary_sizer.Add(count_grid_sizer, 0)

        summary_sizer.Add(wx.StaticText(
            self.page_processing, label="Steps points"), 0)
        self.step_points = wx.StaticText(
            self.page_processing, label='-', size=(200, -1))
        summary_sizer.Add(self.step_points)

        summary_sizer.Add(wx.StaticText(
            self.page_processing, label="Step heights"))
        self.step_heights = wx.StaticText(
            self.page_processing, label='-', size=(200, -1))
        summary_sizer.Add(self.step_heights)

        page_part_menu_sizer.Add(summary_sizer, 0, wx.EXPAND)

        page_part_main_box.Add(page_part_menu_sizer, 0, wx.EXPAND)

        self.page_processing.SetSizer(page_part_main_box, wx.EXPAND)

        defaults = [('method', self.bg_method, 'Pxls sum'),
                    # ('background', self.background, False),
                    ('median_offset', self.offset, int(
                        self.default_offset.GetValue())),
                    ('filter', self.median_filter, False),
                    ('kde', self.kde, False),
                    ('kde_fix', self.kde_fix, False),
                    ('bins', self.bins, int(self.default_bins.GetValue())),
                    ('hist', self.hist, False),
                    ('hist_fix', self.hist_fix, False),
                    ('split_points', self.split_points, ''),
                    ('manual_split', self.manual_split, False)]

        self.particle_widgets = []

        for key, widget, value in defaults:
            gwidget = GWidget(widget, key, value)
            gwidget.setValue(value)
            self.widgets[widget.GetId()] = gwidget
            self.particle_widgets.append(gwidget)

    def generate_photobleaching_page(self):
        self.page_photobleaching = wx.Panel(self.nb)
        self.nb.AddPage(self.page_photobleaching, "Photobleaching")

        self.pb_files_box = wx.StaticBox(
            self.page_photobleaching, 1, " Target Files ", size=(200, -1))
        self.pb_box_sizer = wx.StaticBoxSizer(self.pb_files_box, wx.VERTICAL)
        self.pb_list_spanel = scrolled.ScrolledPanel(self.page_photobleaching)
        self.pb_list_spanel.SetAutoLayout(1)
        self.pb_list_spanel.SetupScrolling()
        self.pb_list_spanel_sizer = wx.BoxSizer(wx.VERTICAL)
        self.pb_list_spanel.SetSizer(self.pb_list_spanel_sizer)
        self.pb_box_sizer.Add(self.pb_list_spanel, 1, wx.EXPAND)

        pb_summary_box = wx.StaticBox(
            self.page_photobleaching, 0, " Photobleaching summary ")
        pb_summary_box_sizer = wx.StaticBoxSizer(pb_summary_box, wx.HORIZONTAL)
        pb_summary_box_sizer.Add(wx.StaticText(
            self.page_photobleaching, label="Number of particles: "), 0)
        self.pbleaching_part_count = wx.StaticText(
            self.page_photobleaching, label="-")
        pb_summary_box_sizer.Add(self.pbleaching_part_count, 1)

        pb_export_btn = wx.Button(
            self.page_photobleaching, label="Export", size=(80, -1))
        pb_export_btn.Bind(
            wx.EVT_BUTTON, self.onExportPhotobleaching, id=pb_export_btn.GetId())

        pb_update_btn = wx.Button(
            self.page_photobleaching, label="Update", size=(80, -1))
        pb_update_btn.Bind(
            wx.EVT_BUTTON, self.onUpdatePhotobleaching, id=pb_update_btn.GetId())

        pb_upbar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pb_upbar_sizer.Add(pb_update_btn, 1, wx.ALIGN_CENTER_VERTICAL)
        pb_upbar_sizer.Add(pb_summary_box_sizer, 1, wx.EXPAND)
        pb_upbar_sizer.Add(pb_export_btn, 1, wx.ALIGN_CENTER_VERTICAL)

        pb_left_box = wx.BoxSizer(wx.VERTICAL)
        pb_left_box.AddSpacer(5)
        pb_left_box.Add(self.pb_box_sizer, 1, wx.EXPAND)

        self.pb_right_sizer = wx.BoxSizer(wx.VERTICAL)
        self.pb_right_sizer.Add(pb_upbar_sizer, 0, wx.EXPAND)
        self.pb_right_sizer.AddSpacer(10)

        self.pbleaching_panel = wx.Panel(self.page_photobleaching)
        self.pb_right_sizer.Add(self.pbleaching_panel, 0)

        tab_pb_main_box = wx.BoxSizer(wx.HORIZONTAL)
        tab_pb_main_box.Add(pb_left_box, 0, wx.EXPAND)
        tab_pb_main_box.Add(self.pb_right_sizer, 1, wx.EXPAND)
        self.page_photobleaching.SetSizer(tab_pb_main_box)

    def generate_brightness_page(self):
        self.page_brightness = wx.Panel(self.nb)
        self.nb.AddPage(self.page_brightness, "Brightness")

        self.ba_files_box = wx.StaticBox(
            self.page_brightness, 1, " Target Files ", size=(200, -1))
        self.ba_box_sizer = wx.StaticBoxSizer(self.ba_files_box, wx.VERTICAL)
        self.ba_list_spanel = scrolled.ScrolledPanel(self.page_brightness)
        self.ba_list_spanel.SetAutoLayout(1)
        self.ba_list_spanel.SetupScrolling()
        self.ba_list_spanel_sizer = wx.BoxSizer(wx.VERTICAL)
        self.ba_list_spanel.SetSizer(self.ba_list_spanel_sizer)
        self.ba_box_sizer.Add(self.ba_list_spanel, 1, wx.EXPAND)

        ba_summary_grid_sizer = wx.FlexGridSizer(4, 2, 2, 2)
        ba_summary_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Number of particles: "), 0)
        self.ba_particles_count = wx.StaticText(
            self.page_brightness, label="-")
        ba_summary_grid_sizer.Add(self.ba_particles_count, 1)
        ba_summary_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Average accuracy: "), 0)
        self.ba_accuracy_avg = wx.StaticText(self.page_brightness, label="-")
        ba_summary_grid_sizer.Add(self.ba_accuracy_avg, 1)
        ba_summary_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Number of monomers: "), 0)
        self.ba_monomers_count = wx.StaticText(self.page_brightness, label="-")
        ba_summary_grid_sizer.Add(self.ba_monomers_count, 1)
        ba_summary_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Monomers intensity: "), 0)
        self.ba_monomer_intensity = wx.StaticText(
            self.page_brightness, label="-")
        ba_summary_grid_sizer.Add(self.ba_monomer_intensity, 1)

        ba_summary_box = wx.StaticBox(self.page_brightness, 0, " Summary ")
        ba_summary_box_sizer = wx.StaticBoxSizer(ba_summary_box, wx.HORIZONTAL)
        ba_summary_box_sizer.Add(ba_summary_grid_sizer, 0, wx.EXPAND)

        ba_export1_btn = wx.Button(
            self.page_brightness, label="Export first frame intensities", size=(80, -1))
        ba_export1_btn.Bind(
            wx.EVT_BUTTON, self.onExportBrightness1, id=ba_export1_btn.GetId())

        ba_export2_btn = wx.Button(
            self.page_brightness, label="Export monomer intensities", size=(80, -1))
        ba_export2_btn.Bind(
            wx.EVT_BUTTON, self.onExportBrightness2, id=ba_export2_btn.GetId())

        ba_export3_btn = wx.Button(
            self.page_brightness, label="Export in separate lines", size=(80, -1))
        ba_export3_btn.Bind(
            wx.EVT_BUTTON, self.onExportBrightness3, id=ba_export3_btn.GetId())

        ba_export_csv = wx.Button(
            self.page_brightness, label="Export as csv", size=(80, -1))
        ba_export_csv.Bind(
            wx.EVT_BUTTON, self.onExportcsv, id=ba_export_csv.GetId())

        ba_results_btn = wx.Button(
            self.page_brightness, label="Show results", size=(80, -1))
        ba_results_btn.Bind(
            wx.EVT_BUTTON, self.onPlotResults, id=ba_results_btn.GetId())

        ba_update_btn = wx.Button(
            self.page_brightness, label="Update", size=(80, -1))
        ba_update_btn.Bind(
            wx.EVT_BUTTON, self.onUpdateBrightness, id=ba_update_btn.GetId())

        opt_box = wx.StaticBox(self.page_brightness, 0, " Options ")
        opt_box_sizer = wx.StaticBoxSizer(opt_box, wx.VERTICAL)
        opt_grid_sizer = wx.FlexGridSizer(4, 2, 2, 2)
        opt_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Label efficiency"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ba_efficiency = wx.TextCtrl(
            self.page_brightness, size=(80, -1), value='0.98')
        opt_grid_sizer.Add(self.ba_efficiency, 0)
        opt_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Brightness:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ba_brightness_method = wx.ComboBox(self.page_brightness, value=self.methods_key['sum'],
                                                choices=[self.methods_key['sum'],
                                                         self.methods_key['2dgintegral'],
                                                         self.methods_key['2dgheight']])
        opt_grid_sizer.Add(self.ba_brightness_method, 0)
        opt_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Background:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ba_bckg_method = wx.ComboBox(self.page_brightness, value=self.methods_key['local'],
                                          choices=[self.methods_key['local'], self.methods_key['global']])
        opt_grid_sizer.Add(self.ba_bckg_method, 0)
        opt_grid_sizer.Add(wx.StaticText(
            self.page_brightness, label="Fitting:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ba_fitting_method = wx.ComboBox(self.page_brightness, value=self.methods_key['gaussfit'],
                                             choices=[self.methods_key['gaussfit'], self.methods_key['pdffit']])
        opt_grid_sizer.Add(self.ba_fitting_method, 0)
        opt_box_sizer.Add(opt_grid_sizer, 0, wx.EXPAND)

        ba_upbar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        ba_upbar_sizer.Add(opt_box_sizer, 1)
        ba_upbar_sizer.Add(ba_update_btn, 1, wx.ALIGN_CENTER_VERTICAL)
        ba_upbar_sizer.Add(ba_summary_box_sizer, 1, wx.EXPAND)
        ba_upbar_sizer.Add(ba_export1_btn, 1, wx.ALIGN_CENTER_VERTICAL)
        ba_upbar_sizer.Add(ba_export2_btn, 1, wx.ALIGN_CENTER_VERTICAL)
        ba_upbar_sizer.Add(ba_export3_btn, 1, wx.ALIGN_CENTER_VERTICAL)
        ba_upbar_sizer.Add(ba_export_csv, 1, wx.ALIGN_CENTER_VERTICAL)

        self.ba_right_sizer = wx.BoxSizer(wx.VERTICAL)
        self.ba_right_sizer.Add(ba_upbar_sizer, 0, wx.EXPAND)
        self.ba_right_sizer.AddSpacer(10)

        self.brightness_panel = wx.Panel(self.page_brightness)
        self.ba_right_sizer.Add(self.brightness_panel, 0)

        ba_left_box = wx.BoxSizer(wx.VERTICAL)
        ba_left_box.AddSpacer(300)
        ba_left_box.Add(self.ba_box_sizer, 1, wx.EXPAND)

        tab_brightness_main_box = wx.BoxSizer(wx.HORIZONTAL)
        tab_brightness_main_box.Add(ba_left_box, 0, wx.EXPAND)
        tab_brightness_main_box.Add(self.ba_right_sizer, 1, wx.EXPAND)
        self.page_brightness.SetSizer(tab_brightness_main_box)

        defaults = [('ba_efficiency', self.ba_efficiency, 0.98),
                    ('ba_method', self.ba_brightness_method,
                     self.methods_key['sum']),
                    ('ba_bckg_method', self.ba_bckg_method, 'Local'),
                    ('ba_fitting_method', self.ba_fitting_method, 'Multiple Gaussians')]

        self.ba_widgets = []

        for key, widget, value in defaults:
            gwidget = GWidget(widget, key, value)
            gwidget.setValue(value)
            self.widgets[widget.GetId()] = gwidget
            self.ba_widgets.append(gwidget)

    def generate_stats_page(self):
        self.page_stats = wx.Panel(self.nb)
        self.nb.AddPage(self.page_stats, "Stats")
        self.tab_stats_main_box = wx.BoxSizer(wx.HORIZONTAL)
        self.stats_panel = wx.Panel(self.page_stats)
        self.stats_panel_sizer = wx.BoxSizer(wx.VERTICAL)
        self.stats_panel_sizer.Add(self.stats_panel, 1, wx.EXPAND)
        self.tab_stats_main_box.Add(self.stats_panel_sizer, 1, wx.EXPAND)
        self.page_stats.SetSizer(self.tab_stats_main_box)

    def generate_analysis_page(self):
        self.page_analysis = wx.Panel(self.nb)
        left_box_analysis = wx.BoxSizer(wx.VERTICAL)
        self.folder_box_loadcsv = wx.StaticBox(
            self.page_analysis, 0, " Load csv ")
        folder_sizer_loadcsv = wx.StaticBoxSizer(self.folder_box_loadcsv, wx.VERTICAL)
        self.folder_path_loadcsv = wx.TextCtrl(self.page_analysis, size=(400, -1))
        folder_sizer_loadcsv.Add(self.folder_path_loadcsv, 0)
        dir_btn_loadcsv = wx.Button(self.page_analysis, label='Browse')
        dir_btn_loadcsv.Bind(wx.EVT_BUTTON, self.onSelectDataFolder_loadcsv,
                     id=dir_btn_loadcsv.GetId())
        folder_sizer_loadcsv.Add(dir_btn_loadcsv, 0, wx.ALIGN_RIGHT)
        left_box_analysis.Add(folder_sizer_loadcsv, 0)
        self.page_analysis.SetSizer(left_box_analysis, wx.EXPAND)
        self.nb.AddPage(self.page_analysis, "Analysis")

        self.csv_files_box = wx.StaticBox(
            self.page_analysis, 1, " Target Files ", size=(400, -1))
        self.folder_sizer_csvfiles = wx.StaticBoxSizer(self.csv_files_box, wx.VERTICAL)
        left_box_analysis.AddSpacer(10)
        left_box_analysis.Add(self.folder_sizer_csvfiles, 0)
        self.csv_list_spanel = scrolled.ScrolledPanel(self.page_analysis)
        self.csv_list_spanel.SetAutoLayout(1)
        self.csv_list_spanel.SetupScrolling()
        self.csv_list_spanel_sizer = wx.BoxSizer(wx.VERTICAL)
        self.csv_list_spanel.SetSizer(self.csv_list_spanel_sizer)
        self.folder_sizer_csvfiles.Add(self.csv_list_spanel, 1, wx.EXPAND)

        analysis_results_btn = wx.Button(
            self.page_analysis, label="Show results", size=(80, -1))
        analysis_results_btn.Bind(
            wx.EVT_BUTTON, self.onPlotResults, id=analysis_results_btn.GetId())
        left_box_analysis.Add(analysis_results_btn, 0)

    def generate_gui_help(self):
        self.folder_box.SetToolTipString(
            "Folder containing all the videos (.tif) ")
        self.default_min_sigma.SetToolTipString("The minimum standard deviation for Gaussian Kernel. "
                                                "Keep this low to detect smaller blobs")
        self.default_max_sigma.SetToolTipString("The maximum standard deviation for Gaussian Kernel. "
                                                "Keep this high to detect larger blobs.")
        self.default_threshold.SetToolTipString("The absolute lower bound for scale space maxima. "
                                                "Local maxima smaller than thresh are ignored. "
                                                "Reduce this to detect blobs with less intensities.")
        self.default_offset.SetToolTipString(
            "Size of the median filter window")
        self.setting_box.SetToolTipString("Checked files will be used as calibration for brightness analysis. "
                                          "Unchecked files will be used for both, photobleaching and brightness analysis")
        self.replicate_combo.SetToolTipString("Video to be processed")
        self.dog_box.SetToolTipString(
            "Automatic detection of bright spots using the Difference of Gaussians method")
        self.min_sigma.SetToolTipString("The minimum standard deviation for Gaussian Kernel. "
                                        "Keep this low to detect smaller blobs")
        self.max_sigma.SetToolTipString("The maximum standard deviation for Gaussian Kernel. "
                                        "Keep this high to detect larger blobs.")
        self.threshold.SetToolTipString("The absolute lower bound for scale space maxima. "
                                        "Local maxima smaller than thresh are ignored. "
                                        "Reduce this to detect blobs with less intensities.")
        self.detect_btn.SetToolTipString(
            "Detect particles using the Difference of Gaussians method")
        self.manual_roi_activation.SetToolTipString(
            "Activate manual manipulation, allows to add or remove particles")
        self.add_opt.SetToolTipString(
            "Add a particle. Click in the image to set coordinates")
        self.add_type.SetToolTipString("Type of the roi to be inserted")
        self.add_radius.SetToolTipString("Radius of the roi to be inserted")
        self.filter_box.SetToolTipString(
            "Filter the signal with a median filter")
        self.offset.SetToolTipString("Size of the median filter window")
        self.kde_box.SetToolTipString("Step detection from peaks in the KDE. ")
        self.hist_box.SetToolTipString(
            "Step detection from peaks in the histogram. ")
        self.pb_files_box.SetToolTipString(" Checked files to be analyzed. ")
        self.ba_files_box.SetToolTipString(" Checked files to be analyzed. ")
        # Todo: complete the help

# Menu events

    def onSave(self, event):
        stch_analysis = {}
        if self.stch_analysis:
            self.update_status_bar("Saving project ...")
            stch_analysis['gui_values'] = self.stch_analysis.gui_values

            replicates = []
            for repl in self.stch_analysis.replicates:
                replicate = {}
                replicate['gui_values'] = repl.gui_values
                replicate['video_name'] = repl.video_name
                replicate['summary'] = repl.summary
                particles = []
                for part in repl.particles:
                    particle = {}
                    particle['x'] = part.x
                    particle['y'] = part.y
                    particle['r'] = part.r
                    particle['localization'] = part.localization
                    particle['gui_values'] = part.gui_values
                    signals = {}
                    for key, signal in part.signals.items():
                        signals[key] = {'frame0_raw_brightness': signal.frame0_raw_brightness,
                                        'frame0_local_bckg': signal.frame0_local_bckg,
                                        'frame0_local_bckg_avg': signal.frame0_local_bckg_avg,
                                        'steps_count': signal.steps_count,
                                        'step_heights': signal.step_heights,
                                        'step_points': signal.step_points}
                    particle['signals'] = signals
                    particles.append(particle)
                replicate['particles'] = particles

                bckg_rois = []
                for bckg in repl.bckg_rois:
                    bckg_roi = {}
                    bckg_roi['r1'] = bckg.r1
                    bckg_roi['r2'] = bckg.r2
                    bckg_roi['c1'] = bckg.c1
                    bckg_roi['c2'] = bckg.c2
                    bckg_rois.append(bckg_roi)
                replicate['bckg_rois'] = bckg_rois

                replicates.append(replicate)
            stch_analysis['replicates'] = replicates

            wildcard = "Project source (*.json)|*.json"
            dlg = wx.FileDialog(self, "Select a project",
                                style=wx.FD_SAVE, wildcard=wildcard)
            if dlg.ShowModal() == wx.ID_OK:
                file = dlg.GetPath()
                relative_path = relpath(self.stch_analysis.folder_path, file)
                stch_analysis['folder_path'] = relative_path
                with open(file, 'w') as outfile:
                    try:
                        json.dump(stch_analysis, outfile)
                    except:
                        wx.MessageBox(
                            'An error ocurred. Project was not saved!', 'Info', wx.OK | wx.ICON_ERROR)
                        self.update_status_bar("Project was not saved!")
                        return
                self.update_status_bar("Project saved!")

    def onLoad(self, event):
        """
        Loading a project
        :param event:
        :return:
        """
        # Loading file dialog
        wildcard = "Project source (*.json)|*.json"
        dlg = wx.FileDialog(self, "Select a project",
                            style=wx.FD_OPEN, wildcard=wildcard)
        if dlg.ShowModal() == wx.ID_OK:
            file_path = dlg.GetPath()
            with open(file_path, 'r') as infile:
                try:
                    stch_analysis = json.load(infile)
                except:
                    wx.MessageBox('Incorrect project file',
                                  'Info', wx.OK | wx.ICON_ERROR)
                    return
            #self.stch_analysis = StchAnalysis(stch_analysis['folder_path'])
            #folder_path = os.path.dirname(file_path)+'/'
            folder_path = realpath(file_path+'/'+stch_analysis['folder_path'])

            analysis_gui_values = {}
            for w in self.analysis_widgets:
                analysis_gui_values[w.key] = w.getValue()

            sequence_gui_values = {}
            for w in self.sequence_widgets:
                sequence_gui_values[w.key] = w.getValue()

            particle_gui_values = {}
            for w in self.particle_widgets:
                particle_gui_values[w.key] = w.getValue()

            self.stch_analysis = StchAnalysis(folder_path,
                                              analysis_gui_values,
                                              sequence_gui_values,
                                              particle_gui_values)
            self.stch_analysis.gui_values = stch_analysis['gui_values']

            # Filling the app objects
            for i, repl in enumerate(stch_analysis['replicates']):
                #  This is for compatibility
                replicate = StchSequence(self.stch_analysis,
                                         repl['video_name'],
                                         repl['gui_values'],
                                         particle_gui_values)
                # This is the one that should be, when only names are saved
                # replicate = StchSequence(folder_path + repl['video_name'])
                fields = ['gui_values', 'video_name', 'summary', ]
                for field in fields:
                    setattr(replicate, field, repl[field])

                for part in repl['particles']:
                    particle = Particle(
                        part['x'], part['y'], part['r'], replicate, part['gui_values'])
                    fields = ['localization']
                    for field in fields:
                        setattr(particle, field, part[field])

                    for key, signal in part['signals'].items():
                        fields = ['frame0_raw_brightness', 'frame0_local_bckg', 'frame0_local_bckg_avg',
                                  'steps_count', 'step_points', 'step_heights']
                        for field in fields:
                            setattr(particle.signals[key],
                                    field, signal[field])
                        particle.signals[key].parent = particle

                    replicate.particles.append(particle)

                for bckg_roi in repl['bckg_rois']:
                    replicate.add_bckg_roi_from_coords(bckg_roi['r1'],
                                                       bckg_roi['r2'],
                                                       bckg_roi['c1'],
                                                       bckg_roi['c2'])

                self.stch_analysis.replicates[i] = replicate

            # Printing projects setting information
            self.folder_path.SetValue(self.stch_analysis.folder_path)
            for gwidget in self.analysis_widgets:
                gwidget.setValue(self.stch_analysis.gui_values[gwidget.key])

            self.update_project_files()

            self.current_sequence = self.stch_analysis.get_replicate_by_index(
                0)
            self.current_particle = None
            self.nb.ChangeSelection(0)

            self.update_status_bar('Project opened!')

    def onClose(self, event):
        # Todo: ask for saving the project
        self.stch_analysis = None
        self.current_sequence = None
        self.current_particle = None
        self.nb.ChangeSelection(0)
        self.folder_path.SetValue("")
        self.update_project_files()
        self.update_status_bar('Project closed')

    def onAbout(self, event):
        description = """
"""

        licence = """Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THE SOFTWARE."""

        info = wx.AboutDialogInfo()

        info.SetIcon(wx.Icon('icon.png', wx.BITMAP_TYPE_PNG))
        info.SetName('Stoichiometry Analysis')
        info.SetVersion(__version__)
        info.SetDescription(description)
        info.SetCopyright(__copyright__)
        info.SetWebSite(
            'http://www.ifib.uni-tuebingen.de/research/garcia-saez.html')
        info.SetLicence(licence)
        info.AddDeveloper(__author__)
        # info.AddDocWriter(__author__)
        # info.AddArtist('')
        # info.AddTranslator('')

        wx.AboutBox(info)


# Common events and methods

    def onChangeTab(self, event):
        self.update_status_bar("Switching tab ...")
        if event.GetSelection() == 0:
            self.update_status_bar("")

        elif event.GetSelection() == 1:  # Moving Image tab
            self.update_image_tab()
            if self.current_sequence:
                self.update_status_bar("{} video loaded!".format(
                    basename(self.current_sequence.video_name)))
            else:
                self.update_status_bar("")

        elif event.GetSelection() == 2:
            self.update_particle_tab()
            self.update_status_bar("")

        elif event.GetSelection() == 3:
            self.update_photobleaching_tab()
            self.update_status_bar("")

        elif event.GetSelection() == 4:
            self.update_brightness_tab()
            self.update_status_bar("")

        elif event.GetSelection() == 5:
            self.update_stats_tab()
            self.update_status_bar("")

    def update_status_bar(self, text):
        try:
            self.status_bar.SetStatusText(text)
        except:
            pass

    def update_gui_value_from_event(self, target, event):
        if target:
            id_key = event.EventObject.GetId()
            gwidget = self.widgets[id_key]
            target.gui_values[gwidget.key] = gwidget.getValue()

    def update_gui_value(self, target, key, value):
        target[key] = value

# Setting s Events and methods

    def onSelectDataFolder(self, event, checked):
        dlg = wx.DirDialog(self, "Select a folder")
        if dlg.ShowModal() == wx.ID_OK:
            dir_path = dlg.GetPath()
            self.folder_path.SetValue(dir_path)

            images = glob.glob(dir_path + '/*.tif')
            if images or checked:
                self.update_status_bar('Loading data folder ...')

                analysis_gui_values = {}
                for w in self.analysis_widgets:
                    analysis_gui_values[w.key] = w.getValue()

                sequence_gui_values = {}
                for w in self.sequence_widgets:
                    sequence_gui_values[w.key] = w.getValue()

                particle_gui_values = {}
                for w in self.particle_widgets:
                    particle_gui_values[w.key] = w.getValue()

                self.update_status_bar('Data folder succefully loaded!')
                if images and not checked:
                    self.stch_analysis = analysis.StchAnalysis(self.folder_path.GetValue(), self.folder_path_cal.GetValue(),
                                                            False,
                                                            analysis_gui_values,
                                                            sequence_gui_values,
                                                            particle_gui_values)
                    self.current_sequence = self.stch_analysis.get_replicate_by_index(
                        0)
                    self.current_particle = None
                    self.update_project_files()
                    self.update_brightness_tab()
            else:
                self.folder_path.SetValue('')
                wx.MessageBox('No TIFF images in selected folder',
                              'Info', wx.OK | wx.ICON_ERROR)
                self.update_status_bar('')

        dlg.Destroy()

    def onSelectDataFolder_loadcsv(self, event):
        dlg = wx.FileDialog(self, "Select file(s)", "", "",
                                   "*.csv", wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            dir_path = dlg.GetPaths()
            print(dir_path)
            self.update_status_bar('Csv files succefully loaded!')

            for cb in self.csv_list_spanel_sizer.GetChildren():  # Removing previous list
                self.csv_list_spanel_sizer.Hide(cb.Window)
                #self.scrolled_panel_sizer.Remove(cb.Window)
                # cb.Window.destroy()
            self.stch_analysis = analysis.StchAnalysis(self)
            for csvfile in dir_path:
                name = csvfile.split("\\")[-1]
                cb = wx.CheckBox(self.csv_list_spanel, label=name)
                cb.Bind(wx.EVT_CHECKBOX, self.onCalibration, id=cb.GetId())
                cb.SetValue(True)
                print(name)
                self.csv_list_spanel_sizer.Add(cb, 0)

                #read in data in dictionary
                with open(csvfile, newline='') as csv_file:
                    reader = csv.reader(csv_file, delimiter=',', quotechar='|')
                    for r,row in enumerate(reader):
                        #calibrations
                        if r == 1:
                            print(row)
                            self.calibrations.append(float(row[0].split(',')[10]))
                        #intensity data
                        if r > 3:
                            for i,key in enumerate(self.fieldnames):
                                try:
                                    print(row)
                                    self.data[key].append(row[0].split(',')[i])
                                except IndexError:
                                    continue
            self.csv_list_spanel_sizer.Layout()
            self.folder_sizer_csvfiles.Layout()
            self.page_analysis.Layout()

    #select folder for calibration files
    def onSelectCalibrationFolder(self, event):
        dlg = wx.DirDialog(self, "Select a folder")
        if dlg.ShowModal() == wx.ID_OK:
            dir_path = dlg.GetPath()
            self.folder_path_cal.SetValue(dir_path)

        dlg.Destroy()

    #select folder to generate tifs out of lifs
    def onSelectDataFolder_ltf(self, event, checked):
        dlg = wx.DirDialog(self, "Select a folder")
        if dlg.ShowModal() == wx.ID_OK:
            dir_path = dlg.GetPath()
            self.folder_path_ltf.SetValue(dir_path)
            
            #generate tif files from available lif files:
            analysis.process_folder(dir_path + "\\", checked)
            print("tif files generates sucessfully !")
        dlg.Destroy()

    #start analysis with paths from data and calibration
    def onStart(self, event):
        rootdir = self.folder_path.GetValue()
        stoichiometries = []    #kinetics of stoichiometry
        n_particles = []    #kinetics of number of particles per image (todo: per mitochondrial area)
        datas = []
        calibrated = False
        for pos,file in enumerate(natsorted(os.listdir(rootdir))):
            d = os.path.join(rootdir, file)
            if os.path.isdir(d):
                print(d)
                istif = glob.glob(d + '/*.tif')
                if istif:
                    self.update_status_bar('Loading data folder ...')

                    analysis_gui_values = {}
                    for w in self.analysis_widgets:
                        analysis_gui_values[w.key] = w.getValue()

                    sequence_gui_values = {}
                    for w in self.sequence_widgets:
                        sequence_gui_values[w.key] = w.getValue()

                    particle_gui_values = {}
                    for w in self.particle_widgets:
                        particle_gui_values[w.key] = w.getValue()

                    self.stch_analysis = analysis.StchAnalysis(d, self.folder_path_cal.GetValue(),
                                                            calibrated,
                                                            analysis_gui_values,
                                                            sequence_gui_values,
                                                            particle_gui_values)
                    self.current_sequence = self.stch_analysis.get_replicate_by_index(
                        0)
                    self.current_particle = None
                    self.update_project_files()

                    self.update_status_bar('Data folder succefully loaded! Starting analysis ...')
                    self.iterate()
                    self.update_brightness_tab()
                    int_data = self.analyzeResults(calibrated)
                    stoichiometry = []
                    n_particle = []
                    for t,dat in enumerate(int_data):
                        stoichiometry.append(np.mean(np.array(dat))/self.stoich_calibration*32)
                        n_particle.append(len(dat))
                        for index,intensity in enumerate(dat):
                            self.data[self.fieldnames[0]].append(rootdir.split("\\")[-1].replace(" ", ""))
                            self.data[self.fieldnames[1]].append(int(pos))
                            self.data[self.fieldnames[2]].append(int(t))
                            self.data[self.fieldnames[3]].append(int(index))
                            self.data[self.fieldnames[4]].append(float(intensity))
                    datas.append(int_data)
                    stoichiometries.append(stoichiometry)
                    n_particles.append(n_particle)
                    calibrated = True
                else:
                    self.folder_path.SetValue('')
                    wx.MessageBox('No TIFF images in selected sub-folder',
                                'Info', wx.OK | wx.ICON_ERROR)
                    self.update_status_bar('')
        
        #we now have an array with a=[[1(t=1),1(t=2),...],[2(t=1),2(t=2),...],...] and we want [[1(t=1),2(t=1),...],[1(t=2),2[t=2],...],...] to be able to perform np.mean(i for i in a)
        s_reshaped = np.array(list(zip(*stoichiometries)))
        s_kinetics = np.array(list(zip(*[(np.mean(i),np.std(i)) for i in s_reshaped])))

        npart_reshaped = np.array(list(zip(*n_particles)))
        npart_kinetics = np.array(list(zip(*[(np.mean(i),np.std(i)) for i in npart_reshaped])))

        data_reshaped = np.array(list(zip(*datas)))

        self.analyzed_data["time"] = np.arange(0, len(s_kinetics[0])*10, 10)
        self.analyzed_data["stoichiometry"] = s_kinetics[0]
        self.analyzed_data["stoichiometry_err"] = s_kinetics[1]
        self.analyzed_data["nparticles"] = npart_kinetics[0]
        self.analyzed_data["nparticles_err"] = npart_kinetics[1]

        print(self.analyzed_data)

        fig, axis = plt.subplots(2,3)
        axis[0][0].set_ylabel("counts")
        axis[0][0].set_xlabel("brightness (mol. units)")
        axis[0][0].set_ylim(0,500)
        axis[0][0].set_xlim(0,600)

        axis[0][1].set_ylabel("stoichiometry")
        axis[0][1].set_xlabel("time after tmre loss (min)")
        axis[0][1].set_xlim(0,110)
        axis[0][1].set_ylim(0,200)

        axis[0][2].set_ylabel("# of spots per area")
        axis[0][2].set_xlabel("time after tmre loss (min)")
        axis[0][2].set_xlim(0,110)
        axis[0][2].set_ylim(0,1000)
        
        dlg = wx.FileDialog(self, "Select a file: ", style=wx.FD_SAVE)

        #plot cumulative distribution of all particles from all samples at time where stoichiometry is max
        #first find maximum
        max_value = max(s_kinetics[0])
        #find index of maximum (if there are multiple take the first one)
        max_index = int([index for index, item in enumerate(s_kinetics[0]) if item == max_value][0])
        print("maximum oligomerization after " + str(max_index*10) + " minutes")
        data_cumulative = []
        for i,p_int in enumerate(datas):
            data_cumulative.append(p_int[max_index]/self.stoich_calibration*32)
        flat_data_cumulative = [num for sublist in data_cumulative for num in sublist]
        axis[0][0].hist(np.array(flat_data_cumulative, dtype = float), alpha = 0.6, bins=fitting.get_bins_number(np.array(flat_data_cumulative, dtype = float)))   #histogram with brightness per spot
        #axis[0][0].legend()

        axis[0][1].plot(np.arange(0, len(s_kinetics[0])*10, 10), s_kinetics[0])
        axis[0][1].fill_between(np.arange(0, len(s_kinetics[0])*10, 10), s_kinetics[0] - s_kinetics[1], s_kinetics[0] + s_kinetics[1], alpha=0.2)
        axis[0][2].plot(np.arange(0, len(npart_kinetics[0])*10, 10), npart_kinetics[0])
        axis[0][2].fill_between(np.arange(0, len(npart_kinetics[0])*10, 10), npart_kinetics[0] - npart_kinetics[1], npart_kinetics[0] + npart_kinetics[1], alpha=0.2)

        fig.show()
    def update_project_files(self):
        """
        Update the list of files of the project
        :return:
        """
        for cb in self.scrolled_panel_sizer.GetChildren():  # Removing previous list
            self.scrolled_panel_sizer.Hide(cb.Window)
            #self.scrolled_panel_sizer.Remove(cb.Window)
            # cb.Window.destroy()
        if self.stch_analysis:
            #all_replicates = self.stch_analysis.replicates + self.stch_analysis.replicates_cal
            for replicate in self.stch_analysis.replicates:
                name = basename(replicate.video_name)
                replicate = self.stch_analysis.get_replicate_by_name(name, False, replicate.summary['calibration'])
                cb = wx.CheckBox(self.scrolled_panel, label=name)
                cb.Bind(wx.EVT_CHECKBOX, self.onCalibration, id=cb.GetId())
                cb.SetValue(replicate.summary['calibration'])
                self.scrolled_panel_sizer.Add(cb, 0)
            self.scrolled_panel_sizer.Layout()
        self.page_settings.Layout()

    def onCalibration(self, event):
        replicate = self.stch_analysis.get_replicate_by_name(
            event.EventObject.Label, load_video=False)
        replicate.summary['calibration'] = event.EventObject.Value

    def onChangeSettingsField(self, event):
        self.update_gui_value_from_event(self.stch_analysis, event)

# Video events and methods

    def update_image_tab(self):
        if self.current_sequence:
            frame = int(self.frame.GetValue())
            img = self.current_sequence.video[frame]
            if (self.current_sequence.particles or self.current_sequence.bckg_rois):
                if self.discard_preview.GetValue():
                    p = float(self.width_threshold.GetValue())
                    idxs = self.current_sequence.get_outliers_indexes(p)
                    img = self.current_sequence.get_image_with_detected_rois(
                        frame, idxs)
                elif self.show_check_box.GetValue():
                    img = self.current_sequence.get_image_with_detected_rois(
                        frame)

            self.update_image_widgets()
            self.update_image_values()
            self.update_main_image(img, basename(
                self.current_sequence.video_name))
            self.manual_roi_activation.SetValue(False)
            self.add_opt.Disable()
            self.remove_opt.Disable()
        else:
            self.update_main_image(None)

    def update_image_widgets(self):
        self.replicate_combo.Clear()
        self.replicate_combo.AppendItems(self.stch_analysis.names)
        self.replicate_combo.SetValue(self.current_sequence.video_name)
        self.frame.SetRange(0, len(self.current_sequence.video)-1)
        # Check for default_values
        if len(self.current_sequence.particles) == 0:
            self.current_sequence.gui_values['min_sigma'] = self.stch_analysis.gui_values['def_min_sigma']
            self.current_sequence.gui_values['max_sigma'] = self.stch_analysis.gui_values['def_max_sigma']
            self.current_sequence.gui_values['threshold'] = self.stch_analysis.gui_values['def_threshold']
            self.current_sequence.gui_values['roi_radius'] = self.stch_analysis.gui_values['def_roi_radius']
        for w in self.sequence_widgets:
            w.setValue(self.current_sequence.gui_values[w.key])

    def update_image_values(self):
        coord_x = '-'
        coord_y = '-'
        sum_bgth = '-'
        sum_bckg = '-'
        total = 0
        global_bckg_avg = '-'
        local_bckg_avg = '-'
        bg_method = self.methods_val[self.bg_method.GetValue()]
        part_accr_local = '-'
        part_accr_global = '-'
        accr_local_avg = '-'
        accr_global_avg = '-'
        fitted_width = '-'
        fitted_width_avg = '-'
        fitted_width_std = '-'
        if self.current_sequence:
            frame = int(self.frame.GetValue())
            self.current_sequence.update_global_background(frame)
            self.current_sequence.update_summary_values(bg_method)
            if self.current_particle:
                coord_x = str(self.current_particle.x)
                coord_y = str(self.current_particle.y)
                sum_bgth = "{:.2f}".format(
                    self.current_particle.signals[bg_method].frame0_raw_brightness)
                sum_bckg = "{:.2f}".format(
                    self.current_particle.signals[bg_method].frame0_local_bckg_avg)
                if self.current_particle.localization != {}:
                    part_accr_global = "{:.2f}".format(
                        self.current_particle.localization['accuracy_global'])
                    part_accr_local = "{:.2f}".format(
                        self.current_particle.localization['accuracy_local'])
                    fitted_width = "{:.2f}".format(
                        self.current_particle.localization['d'])
            if self.current_sequence.particles:
                total = len(self.current_sequence.particles)
            global_bckg_avg = "{:.2f}".format(
                self.current_sequence.summary['global_bckg_avg'])
            local_bckg_avg = "{:.2f}".format(
                self.current_sequence.summary['local_bckg_avg'])
            accr_local_avg = "{:.2f}".format(
                self.current_sequence.summary['accuracy_local'])
            accr_global_avg = "{:.2f}".format(
                self.current_sequence.summary['accuracy_global'])
            fitted_width_avg = "{:.2f}".format(
                self.current_sequence.summary['fitted_width_avg'])
            fitted_width_std = "{:.2f}".format(
                self.current_sequence.summary['fitted_width_std'])

        self.calibration.SetLabel(
            str(self.current_sequence.summary['calibration']))
        self.coords_x_part.SetLabel(coord_x)
        self.coords_y_part.SetLabel(coord_y)
        self.part_brightness_sum.SetLabel(sum_bgth)
        self.part_brightness_sum_bckg.SetLabel(sum_bckg)
        self.summary_box.SetLabel(" Summary ({})".format(total))
        self.global_bckg_avg.SetLabel(global_bckg_avg)
        self.local_bckg_avg.SetLabel(local_bckg_avg)
        self.part_accr_global.SetLabel(part_accr_global)
        self.part_accr_local.SetLabel(part_accr_local)
        self.accr_global_avg.SetLabel(accr_global_avg)
        self.accr_local_avg.SetLabel(accr_local_avg)
        self.fitted_width.SetLabel(fitted_width)
        self.fitted_width_avg.SetLabel(fitted_width_avg)
        self.fitted_width_std.SetLabel(fitted_width_std)
        self.page_detection.Fit()

    def update_main_image(self, image, title=""):
        self.page_detection_image_sizer.Hide(self.page_detection_image_panel)
       # self.page_detection_image_sizer.Remove(self.page_detection_image_panel)
        self.page_detection_image_panel.Destroy()
        self.page_detection_image_panel = wx.Panel(self.page_detection)

        if image is not None:

            figure = Figure()
            axes = figure.add_subplot(111)
            axes.set_title(title)
            axes.imshow(image, cmap=plt.cm.gray, interpolation="nearest")
            canvas = FigureCanvas(self.page_detection_image_panel, 1, figure)

            def onClick(e):
                x = e.ydata  # coords invertedSetStatusText
                y = e.xdata
                particle_id = self.current_sequence.get_particle_from_coordinates(
                    x, y)
                bckg_idx = self.current_sequence.get_background_from_coordinates(
                    x, y)
                frame = int(self.frame.GetValue())
                if self.manual_roi_activation.GetValue():
                    if self.remove_opt.GetValue():
                        if particle_id is not None:
                            self.current_sequence.particles.pop(particle_id)
                            img = self.current_sequence.get_image_with_detected_rois(
                                frame)
                            self.update_main_image(img)
                            self.current_particle = None
                            self.update_status_bar(
                                " Particle {} removed".format(particle_id))
                        elif bckg_idx is not None:
                            self.current_sequence.bckg_rois.pop(bckg_idx)
                            img = self.current_sequence.get_image_with_detected_rois(
                                frame)
                            self.update_main_image(img)
                            self.update_status_bar(" Background ROI removed")
                    if self.add_opt.GetValue():
                        x = np.round(x).astype(int)  # coords inverted
                        y = np.round(y).astype(int)
                        r = self.add_radius.GetValue()
                        if self.add_type.GetValue() == 'Particle':
                            particle_id = self.current_sequence.add_particle(
                                x, y, r)
                            if particle_id is not None:
                                img = self.current_sequence.get_image_with_detected_rois(
                                    frame)
                                self.current_particle = self.current_sequence.particles[-1]
                                self.show_check_box.SetValue(True)
                                self.update_image_values()
                                self.update_main_image(img)
                                self.update_status_bar(
                                    " Particle {} added at ({}, {})".format(particle_id, y, x))
                            else:
                                self.update_status_bar(" Particle not added ")
                        elif self.add_type.GetValue() == 'Background':
                            res = self.current_sequence.add_bckg_roi_from_width(
                                x, y, r)
                            if res is not None:
                                self.show_check_box.SetValue(True)
                                img = self.current_sequence.get_image_with_detected_rois(
                                    frame)
                                self.update_image_values()
                                self.update_main_image(img)
                                self.update_status_bar(" Background ROI added")
                            else:
                                self.update_status_bar(
                                    " Background ROI not added")
                else:
                    if particle_id is not None:
                        self.current_particle = self.current_sequence.particles[particle_id]
                        self.update_image_values()
                        self.particle_id.SetValue(particle_id)
                    try:
                        self.update_status_bar(
                            " Coordinates: ({:.2f}, {:.2f})".format(y, x))
                    except:
                        pass

            canvas.mpl_connect('button_press_event', onClick)
            try:
                toolbar = NavigationToolbar2Wx(canvas)
                toolbar.Realize()
            except wx._core.wxAssertionError:
                print("canvas error")
            panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
            panel_sizer.Add(canvas, 1, wx.EXPAND)
            panel_sizer.Add(toolbar, 0)
        else:
            panel_sizer = wx.BoxSizer(wx.VERTICAL)
            panel_sizer.Add(wx.Panel(self.page_detection), 1, wx.EXPAND)

        self.page_detection_image_panel.SetSizer(panel_sizer, wx.EXPAND)
        self.page_detection_image_sizer.Add(
            self.page_detection_image_panel, 1, wx.EXPAND)
        self.page_detection.Layout()

    def onChangeImage(self, event):
        """
        What to do when selecting a different image
        :param event:
        :return:
        """
        image_name = self.replicate_combo.GetValue()
        self.update_status_bar("Loading {} video ...".format(image_name))
        self.current_sequence = self.stch_analysis.get_replicate_by_name(
            image_name)
        self.current_particle = None
        self.update_image_tab()
        self.update_status_bar("{} loaded!".format(image_name))

    def onChangeFrame(self, event):
        self.update_gui_value_from_event(self.current_sequence, event)
        self.update_image_tab()

    def onChangeDetectionField(self, event):
        self.update_gui_value_from_event(self.current_sequence, event)

    def onBackground(self, event):
        self.update_status_bar("Generating background ROIs")
        frame = int(self.frame.GetValue())
        width = int(self.add_radius.GetValue())*2
        self.current_sequence.generate_background(frame, width)
        self.show_check_box.SetValue(True)
        self.update_image_tab()
        self.update_status_bar("Background ROIs generated!")

    def onDetect(self, event):
        self.update_status_bar('Detecting particles ...')
        self.current_sequence.particles = []
        res = self.current_sequence.detect_particles(
            int(self.frame.GetValue()))
        self.show_check_box.SetValue(True)
        if self.current_sequence.particles:
            self.current_particle = self.current_sequence.particles[0]
        self.update_image_tab()
        self.update_status_bar("{0} particles detected".format(res))

    def onLocalize(self, event):
        self.update_status_bar('Localizing particles ...')
        frame = int(self.frame.GetValue())
        fit_method = self.fit_method.GetValue()
        self.current_sequence.localize_particles(frame, fit_method)
        self.show_check_box.SetValue(True)
        self.update_image_tab()
        self.update_status_bar("Localization is done!")

    def onIterate(self, event):
        iterate(self)
    
    def iterate(self):
        self.update_status_bar('Iterating ...')
        #all_replicates = self.stch_analysis.replicates + self.stch_analysis.replicates_cal
        for x in self.stch_analysis.replicates:
            yt = basename(x.video_name)

            self.update_status_bar("Loading {} video ...".format(yt))
            self.current_sequence = self.stch_analysis.get_replicate_by_name(
                yt, calibration = x.summary['calibration'])
            self.current_particle = None
            self.update_image_tab()
            self.update_status_bar("{} loaded!".format(yt))

            self.update_status_bar('Detecting particles ...')
            self.current_sequence.particles = []
            res = self.current_sequence.detect_particles(
                int(self.frame.GetValue()))
            self.show_check_box.SetValue(True)
            if self.current_sequence.particles:
                self.current_particle = self.current_sequence.particles[0]
            self.update_image_tab()
            self.update_status_bar("{0} particles detected".format(res))
            print("detected")
            #self.update_status_bar("Generating background ROIs")
            #frame = int(self.frame.GetValue())
            #width = int(self.add_radius.GetValue())*2
            #self.current_sequence.generate_background(frame, width)
            # self.show_check_box.SetValue(True)
            # self.update_image_tab()
            #self.update_status_bar("Background ROIs generated!")
            
            for itnum in range(int(self.iterations.GetValue())):
                try:
                    print(itnum)
                    self.update_status_bar('Iterating ...')
                    frame = int(self.frame.GetValue())
                    fit_method = self.fit_method.GetValue()
                    self.current_sequence.localize_particles(frame, fit_method)
                    self.show_check_box.SetValue(True)
                    self.update_image_tab()

                    #self.update_status_bar('Discarding particles ...')
                    # elf.current_sequence.discard_close_particles()
                    # self.show_check_box.SetValue(True)
                    # self.update_image_tab()
                    #self.update_status_bar("Discarding is done")
                    
                    self.update_status_bar('Discarding particles (Threshold)...')
                    threshold = float(self.width_threshold.GetValue())
                    self.current_sequence.discard_wide_particles(threshold)
                    self.discard_preview.SetValue(False)
                    self.show_check_box.SetValue(True)
                    self.update_image_tab()
                    self.update_status_bar("Discarding is done")

                    self.update_status_bar('Discarding particles (Stack)...')
                    # Gathering the values for all particles
                    to_remove = []
                    for i, particle in enumerate(self.current_sequence.particles):
                        #print("particle %i" %i)
                        particle.signals['sum'].generate_sequence()
                        # Removing condition
                        # Question: what if the maximum is in both, the middle and one of the extremes
                        # discarding first frame
                        frame_values = particle.signals['sum'].base_values[1:]
                        max_value = max(frame_values)
                        if frame_values[0] == max_value or frame_values[-1] == max_value:
                            to_remove.append(i)
                    self.current_sequence.discard_particles(to_remove)
                    self.update_image_tab()
                    self.update_status_bar("Discarding is done...")
                except ValueError:
                    print(str(self.current_sequence.video_name) + " is probably not suitable for analysis because no spots could be detected")
                    continue

            self.update_status_bar("Iteration is done!")

    def onDiscard(self, event):
        self.update_status_bar('Discarding particles ...')
        self.current_sequence.discard_close_particles()
        self.show_check_box.SetValue(True)
        self.update_image_tab()
        self.update_status_bar("Discarding is done")

    def onWDiscard(self, event):
        self.update_status_bar('Discarding particles ...')
        threshold = float(self.width_threshold.GetValue())
        self.current_sequence.discard_wide_particles(threshold)
        self.discard_preview.SetValue(False)
        self.show_check_box.SetValue(True)
        self.update_image_tab()
        self.update_status_bar("Discarding is done")

    def onDiscardPreview(self, event):
        self.update_image_tab()

    def onStackDiscard(self, event):
        self.update_status_bar('Discarding particles ...')
        # Gathering the values for all particles
        to_remove = []
        for i, particle in enumerate(self.current_sequence.particles):
            particle.signals['sum'].generate_sequence()
            # Removing condition
            # Question: what if the maximum is in both, the middle and one of the extremes
            # discarding first frame
            frame_values = particle.signals['sum'].base_values[1:]
            max_value = max(frame_values)
            if frame_values[0] == max_value or frame_values[-1] == max_value:
                to_remove.append(i)
        self.current_sequence.discard_particles(to_remove)
        self.update_image_tab()
        self.update_status_bar("Discarding is done")

    def onManualParticleHandling(self, event):
        if self.manual_roi_activation.GetValue():
            action = "Enable"
        else:
            action = "Disable"
        for w in self.manual_roi_controls:
            getattr(w, action)()

    def onShowParticles(self, event):
        self.update_image_tab()

    def onVisualize(self, event):
        if self.current_sequence and self.current_particle:
            x, y, r = self.current_particle.x, self.current_particle.y, self.current_particle.r
            X, Y = np.meshgrid(np.arange(x-r, x+r+1),
                               np.arange(y-r, y+r+1),
                               indexing='ij')
            frame = int(self.frame.GetValue())
            roi = self.current_sequence.video[frame][X, Y]
            fig = plt.figure()
            fit = self.current_particle.get_fitted(frame)
            if fit is not None:
                ax = fig.add_subplot(1, 2, 1, projection='3d')
                surf = ax.plot_surface(X, Y, roi, rstride=1, cstride=1, cmap=cm.coolwarm,
                                       linewidth=0, antialiased=False)
                plt.title("Raw data: ({}, {})".format(x, y))
                ax = fig.add_subplot(1, 2, 2, projection='3d')
                surf = ax.plot_surface(X, Y, fit, rstride=1, cstride=1, cmap=cm.coolwarm,
                                       linewidth=0, antialiased=False)
                xc = self.current_particle.localization['xc']
                yc = self.current_particle.localization['yc']
                plt.title("2D Gaussian: ({:.2f}, {:.2f})".format(xc, yc))

            else:
                ax = fig.add_subplot(1, 1, 1, projection='3d')
                surf = ax.plot_surface(X, Y, roi, rstride=1, cstride=1, cmap=cm.coolwarm,
                                       linewidth=0, antialiased=False)
                plt.title("Raw data: ({}, {})".format(x, y))

            plt.show()

# Particle events and methods

    def onChangeParticle(self, event):
        self.current_particle = self.current_sequence.particles[self.particle_id.GetValue(
        )]
        self.update_particle_tab()

    def onDeleteParticle(self, event):
        if self.current_particle:
            idx = self.particle_id.GetValue()
            self.current_sequence.particles.pop(idx)
            if idx >= len(self.current_sequence.particles):
                idx -= 1
            if idx >= 0:
                self.current_particle = self.current_sequence.particles[idx]
            else:
                self.current_particle = None
            self.update_particle_tab()

    def onChangeParticleField(self, event):
        self.update_gui_value_from_event(self.current_particle, event)
        self.update_particle_tab()

    def onAutoBin(self, event):
        """
        Compute an apropiate number of bins based on Freedman-Diaconis rule
        :param event:
        :return:
        """
        ba_method = self.methods_val[self.ba_brightness_method.GetValue()]
        if self.current_particle:
            self.current_particle.signals[ba_method].generate_sequence()
            bins = get_bins_number(
                self.current_particle.signals[ba_method].base_values)
            bin_widget = self.widgets[self.bins.GetId()]
            bin_widget.setValue(bins)
            self.current_particle.gui_values[bin_widget.key] = bins
        self.update_particle_tab()

    def onChangeSplitPoints(self, event):
        if self.current_particle:
            try:
                ast.literal_eval('[' + self.split_points.GetValue() + ']')
                split_points_widget = self.widgets[self.split_points.GetId()]
                self.current_particle.gui_values[split_points_widget.key] = self.split_points.GetValue(
                )
            except:
                # Todo: show a message
                pass

    def update_particle_tab(self, trigger=None):
        if self.current_sequence and self.current_sequence.particles:
            if not self.current_particle:
                self.current_particle = self.current_sequence.particles[0]
            self.update_particle_widgets()
            self.update_particle_image(trigger=trigger)
            self.update_particle_values()

    def update_particle_widgets(self):
        self.particle_id.SetRange(0, len(self.current_sequence.particles)-1)
        if not self.current_particle.gui_values['filter']:
            self.current_particle.gui_values['median_offset'] = self.stch_analysis.gui_values['def_median_offset']
        if not (self.current_particle.gui_values['hist'] or self.current_particle.gui_values['hist_fix']):
            self.current_particle.gui_values['bins'] = self.stch_analysis.gui_values['def_bins']
        for w in self.particle_widgets:
            w.setValue(self.current_particle.gui_values[w.key])

    def update_particle_values(self):
        if self.current_particle:
            coord_x = str(self.current_particle.x)
            coord_y = str(self.current_particle.y)
            bg_method = self.methods_val[self.bg_method.GetValue()]
            kde_steps_count = str(
                self.current_particle.signals[bg_method].kde_steps)
            kde_fix_steps_count = str(
                self.current_particle.signals[bg_method].kde_fixed_steps)
            hist_steps_count = str(
                self.current_particle.signals[bg_method].hist_steps)
            hist_fix_steps_count = str(
                self.current_particle.signals[bg_method].hist_fixed_steps)
            manual_steps_count = str(
                self.current_particle.signals[bg_method].manual_steps)
            steps_count = str(
                self.current_particle.signals[bg_method].steps_count)
            step_points = str(
                self.current_particle.signals[bg_method].step_points)
            photon_coef = self.widgets[self.photon_coef.GetId()].getValue()
            step_heights = str(
                [int(x*photon_coef) for x in self.current_particle.signals[bg_method].step_heights])
            if self.current_particle.localization != {}:
                width = "{:.2f}".format(
                    self.current_particle.localization['d'])
            else:
                width = '-'
        else:
            coord_x = '-'
            coord_y = '-'
            width = '-'
            kde_steps_count = '-'
            kde_fix_steps_count = '-'
            hist_steps_count = '-'
            hist_fix_steps_count = '-'
            manual_steps_count = '-'
            steps_count = '-'
            step_points = '-'
            step_heights = '-'

        self.coords_x.SetLabel(coord_x)
        self.coords_y.SetLabel(coord_y)
        self.width.SetLabel(width)
        self.kde_steps_count.SetLabel(kde_steps_count)
        self.kde_fix_steps_count.SetLabel(kde_fix_steps_count)
        self.hist_steps_count.SetLabel(hist_steps_count)
        self.hist_fix_steps_count.SetLabel(hist_fix_steps_count)
        self.manual_steps_count.SetLabel(manual_steps_count)
        self.steps_count.SetLabel(steps_count)
        self.step_points.SetLabel(step_points)
        self.step_points.Wrap(195)
        self.step_heights.SetLabel(step_heights)
        self.step_heights.Wrap(195)
        self.page_processing.Fit()

    def update_particle_image(self, title="", trigger=None):
        """
        Main method for plotting intensities over time
        :param title:
        :param trigger: which checkbox triggered the method
        :return:
        """
        self.page_part_image_sizer.Hide(self.page_part_image_panel)
        # self.page_part_image_sizer.RemoveChild(self.page_part_image_panel)
        self.page_part_image_panel = wx.Panel(self.page_processing)

        if self.current_particle:

            figure = Figure()
            axes = figure.add_axes([0.07, 0.07, 0.9, 0.9])

            bg_method = self.methods_val[self.bg_method.GetValue()]
            gsignal = self.current_particle.signals[bg_method]
            gsignal.generate_signals()

            # (ys, xs, alpha, color)
            signals = [(gsignal.raw_values, 0.3, 'k'),
                       (gsignal.filtered_values, 0.5, 'g'),
                       (gsignal.kde_values, 0.7, 'b'),
                       (gsignal.kde_fixed_values, 1.0, 'b'),
                       (gsignal.hist_values, 0.7, 'r'),
                       (gsignal.hist_fixed_values, 1.0, 'r'),
                       (gsignal.manual_values, 1.0, 'c')]

            for signal in signals:
                if signal[0].any():
                    ys = signal[0]
                    xs = range(len(gsignal.raw_values))
                    alpha = signal[1]
                    color = signal[2]
                    axes.plot(np.multiply(xs, float(self.stch_analysis.gui_values['frame_rate'])),
                              np.multiply(
                                  ys, float(self.stch_analysis.gui_values['photon_coef'])),
                              alpha=alpha, color=color)
                    axes.set_xlabel("Time (miliseconds)")
                    axes.set_ylabel("Photon count")

            if gsignal.kde_values.any() or gsignal.kde_fixed_values.any() or \
                    gsignal.hist_values.any() or gsignal.hist_fixed_values.any():
                inset_plot = figure.add_axes([0.7, 0.7, 0.20, 0.20])
                if gsignal.kde_values.any() or gsignal.kde_fixed_values.any():
                    inset_plot.plot(gsignal.kde_xs,
                                    gsignal.kde_ys,
                                    color='b')
                    for peak in gsignal.kde_peak_values:
                        inset_plot.axvline(peak, ls=':', color='b')
                if gsignal.hist_values.any() or gsignal.hist_fixed_values.any():
                    width = 1.0 * \
                        (gsignal.hist_bin_edges[1] - gsignal.hist_bin_edges[0])
                    center = (
                        gsignal.hist_bin_edges[:-1] + gsignal.hist_bin_edges[1:]) / 2
                    inset_plot.bar(center, gsignal.hist_bin_values,
                                   align='center', width=width, color='r', alpha=0.5)
                    for peak in gsignal.hist_peaks_values:
                        inset_plot.axvline(peak, ls=':', color='r')

            canvas = FigureCanvas(self.page_part_image_panel, 1, figure)

            def onClick(e):
                x = e.ydata  # coords inverted
                y = e.xdata
                try:
                    fr = float(self.stch_analysis.gui_values['frame_rate'])
                    pc = float(self.stch_analysis.gui_values['photon_coef'])
                    self.update_status_bar(
                        " Coordinates: ({:.2f}, {:.2f}), Raw values: ({:.2f}, {:.2f})".format(y, x, y/fr, x/pc))
                except:
                    pass

            canvas.mpl_connect('button_press_event', onClick)

            toolbar = NavigationToolbar2Wx(canvas)
            toolbar.Realize()

            panel_sizer = wx.BoxSizer(wx.VERTICAL)
            panel_sizer.Add(canvas, 1, wx.EXPAND)
            panel_sizer.Add(toolbar, 0)
        else:
            panel_sizer = wx.BoxSizer(wx.VERTICAL)
            panel_sizer.Add(wx.Panel(self.page_processing), 1, wx.EXPAND)

        self.page_part_image_panel.SetSizer(panel_sizer)
        self.page_part_image_sizer.Add(
            self.page_part_image_panel, 1, wx.EXPAND)
        self.page_processing.Layout()

# Photobleaching events and methods

    def onExportPhotobleaching(self, event):
        self.update_status_bar("Exporting photobleaching results ...")
        if self.stch_analysis:
            str = self.stch_analysis.export_photobleaching()
            dlg = wx.FileDialog(self, "Select a file: ", style=wx.FD_SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                file = dlg.GetPath()
                with open(file, 'w') as outfile:
                    outfile.write(str)
                    outfile.close()
        self.update_status_bar("Photobleaching results succefully exported!")

    def onUpdatePhotobleaching(self, event):
        self.update_photobleaching_tab()

    def update_photobleaching_tab(self):
        self.update_file_list(self.page_photobleaching,
                              self.pb_list_spanel,
                              self.pb_list_spanel_sizer,
                              'pb_analysis',
                              self.onPhotobleachingAnalysis,
                              'no_calibration')
        self.update_photobleaching_image()

    def update_photobleaching_image(self):
        self.pb_right_sizer.Hide(self.pbleaching_panel)
        # self.pb_right_sizer.RemoveChild(self.pbleaching_panel)
        self.pbleaching_panel = wx.Panel(self.page_photobleaching)
        panel_sizer = wx.BoxSizer(wx.VERTICAL)

        if self.stch_analysis is not None:
            method = self.methods_val[self.ba_brightness_method.GetValue()]
            self.stch_analysis.photobleaching_analysis(method)
            if self.stch_analysis.pb_count and self.stch_analysis.pb_steps_count:
                figure = Figure()

                # Species counting
                axes = figure.add_subplot(121)

                count = 0
                ys = np.zeros(7)
                for k, v in self.stch_analysis.pb_count.items():
                    if 1 <= k <= 6:
                        ys[k - 1] = v
                        count += v
                    elif k > 6:
                        ys[-1] = v
                        count += v

                self.pbleaching_part_count.SetLabel(str(count))
                xs = np.arange(1, 8)
                axes.set_xticks(xs)
                width = 0.8
                xs = xs-0.4
                axes.bar(xs, ys, width=width, alpha=0.5)
                axes.set_xlim([0.5, 7.5])
                axes.set_title('Photobleaching Analysis')
                axes.set_ylabel('Count')
                axes.set_xlabel('Species (number of steps)')

                # Heights distribution
                axes = figure.add_subplot(122)

                data = [[] for i in range(7)]
                photon_count = self.stch_analysis.gui_values['photon_coef']
                for (specie, step), step_height in self.stch_analysis.pb_steps_count.items():
                    if 1 <= specie <= 6:
                        data[specie -
                             1].extend(np.array(step_height)*photon_count)
                    elif specie > 6:
                        data[-1].extend(np.array(step_height)*photon_count)

                bp = axes.boxplot(data, sym='.', showfliers=False)

                for x, values in enumerate(data):
                    for y in values:
                        axes.plot(x+1, y, 'k.', alpha=0.5)

                axes.set_title('Step heights')
                axes.set_ylabel('Intensity')
                axes.set_xlabel('Species (number of steps)')
                figure.subplots_adjust(
                    left=0.05, right=0.95, top=0.95, bottom=0.1)

                canvas = FigureCanvas(self.pbleaching_panel, 1, figure)
                toolbar = NavigationToolbar2Wx(canvas)
                toolbar.Realize()

                panel_sizer.Add(canvas, 1, wx.EXPAND)
                panel_sizer.Add(toolbar, 0)

            else:
                self.pbleaching_part_count.SetLabel("-")
        else:
            self.pbleaching_part_count.SetLabel("-")

        self.pbleaching_panel.SetSizer(panel_sizer)
        self.pb_right_sizer.Add(self.pbleaching_panel, 1, wx.EXPAND)
        self.page_photobleaching.Layout()

    def onPhotobleachingAnalysis(self, event):
        replicate = self.stch_analysis.get_replicate_by_name(
            event.EventObject.Label, load_video=False)
        replicate.summary['pb_analysis'] = event.EventObject.Value

    def update_file_list(self, page, spanel, spanel_sizer, field, method, filter='all'):
        """
        Update a scrolled panel with files list
        :param page:
        :param scrolled_panel:
        :param field: meaning of the checkbox (calibration, ba_analysis, pb_analysis)
        :param filter: 'all', 'calibration', 'no_calibration'
        :return:
        """
        for cb in spanel_sizer.GetChildren():  # Removing previous list
            spanel_sizer.Hide(cb.Window)
            #spanel_sizer.Remove(cb.Window)
        if self.stch_analysis:
            for replicate in self.stch_analysis.replicates:
                if (filter == 'calibration' and replicate.summary['calibration']) or \
                        (filter == 'no_calibration' and not replicate.summary['calibration']) or \
                        (filter == 'all'):
                    name = basename(replicate.video_name)
                    replicate = self.stch_analysis.get_replicate_by_name(
                        name, load_video=False)
                    cb = wx.CheckBox(spanel, label=name)
                    cb.Bind(wx.EVT_CHECKBOX, method, id=cb.GetId())
                    cb.SetValue(replicate.summary[field])
                    spanel_sizer.Add(cb, 0)
            spanel_sizer.Layout()
        page.Layout()

# Brightness events and methods

    def onExportBrightness1(self, event):
        self.update_status_bar("Exporting brightness analysis results ...")
        str = self.stch_analysis.export_brightness_1()
        dlg = wx.FileDialog(self, "Select a file: ", style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            file = dlg.GetPath()
            with open(file, 'w') as outfile:
                outfile.write(str)
                outfile.close()
        self.update_status_bar(
            "Brightness analysis results succefully exported!")

    def onExportBrightness2(self, event):
        self.update_status_bar("Exporting brightness analysis results ...")
        str = self.stch_analysis.export_brightness_2()
        dlg = wx.FileDialog(self, "Select a file: ", style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            file = dlg.GetPath()
            with open(file, 'w') as outfile:
                outfile.write(str)
                outfile.close()
        self.update_status_bar(
            "Brightness analysis results succefully exported!")

        # export in separate lines
    def onExportBrightness3(self, event):
        self.update_status_bar("Exporting brightness analysis results ...")
        str = self.stch_analysis.export_brightness_3()
        dlg = wx.FileDialog(self, "Select a file: ", style=wx.FD_SAVE)
        # separate into lines:
        splitstr = str.split(",")
        combstr = ""
        for i,replicate in enumerate(self.stch_analysis.replicates):
            for j in range(len(replicate.particles)-1):
                combstr += splitstr[0]
                splitstr.pop(0)
                combstr += " "
            combstr += "\n"
        if dlg.ShowModal() == wx.ID_OK:
            file = dlg.GetPath()
            with open(file, 'w') as outfile:
                outfile.write(combstr)
                outfile.close()
        self.update_status_bar(
            "Brightness analysis results succefully exported!")

    # export dictionary as .csv
    def onExportcsv(self, *params):    #params: samples = [], position (1-10)
        self.update_status_bar("Exporting analysis results ...")
        dlg = wx.FileDialog(self, "Select a file: ", style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            #export raw data
            file = dlg.GetPath()
            with open(file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter = ",")
                options_keys = [*self.stch_analysis.gui_values.keys(), 'calibration']
                writer.writerow(options_keys)
                options_vals = [*self.stch_analysis.gui_values.values(),self.stoich_calibration]
                writer.writerow(options_vals)
                writer.writerow('') #add empty line
                writer.writerow(self.data.keys())
                writer.writerows(zip(*self.data.values()))

            print(self.analyzed_data)
            print(zip(*self.analyzed_data.values()))
            #export analyzed data
            with open(file.split(".")[0] + "_analyzed.csv", 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter = ",")
                writer.writerow(self.analyzed_data.keys())
                writer.writerows(zip(*self.analyzed_data.values()))

        self.update_status_bar(
            "Brightness analysis results succefully exported as csv!")

    def onPlotResults(self, event):
        self.plotResults()

    def plotResults(self):
        fig, axis = plt.subplots(2,3)
        axis[0][0].set_ylabel("counts")
        axis[0][0].set_xlabel("brightness (mol. units)")
        axis[0][0].set_ylim(0,500)
        axis[0][0].set_xlim(0,600)

        axis[0][1].set_ylabel("stoichiometry")
        axis[0][1].set_xlabel("time after tmre loss (min)")
        axis[0][1].set_xlim(0,110)
        axis[0][1].set_ylim(0,200)

        axis[0][2].set_ylabel("# of spots per area")
        axis[0][2].set_xlabel("time after tmre loss (min)")
        axis[0][2].set_xlim(0,110)
        axis[0][2].set_ylim(0,1000)

        axis[1][1].set_ylabel("stoichiometry")
        axis[1][1].set_xlabel("time after tmre loss (min)")
        axis[1][1].set_xlim(0,110)
        axis[1][1].set_ylim(0,200)

        axis[0][0].hist(np.array(self.data['intensity'], dtype = np.float)/self.calibrations[0]*32)

        n_times = int(max(self.data['time']))
        stoichiometry = np.zeros(n_times+1)
        print(n_times)
        #now plot different experiments in different lines
        #1: find all different experiments in dictionary
        experiment_names = []
        for e in self.data['experiment']:
            if not e in experiment_names:
                experiment_names.append(e)
        print(experiment_names)
        #construct nested list
        intensity_experiments = []
        for e in experiment_names:
            intensity_times = []
            for t in range(n_times+1):
                intensity_times.append([])
            #fill lists with intensity values timewise
            for i in range(len(self.data['time'])):
                if self.data["experiment"][i] == e:
                    intensity_times[int(self.data['time'][i])].append(float(self.data['intensity'][i]))
            stoichiometry = np.array(list(zip(*[(np.mean(i),np.std(i)) for i in intensity_times])))/self.calibrations[0]*32
            intensity_experiments.append(intensity_times)
            axis[0][1].plot(np.arange(0, len(stoichiometry[0])*10, 10), stoichiometry[0])
            axis[0][1].fill_between(np.arange(0, len(stoichiometry[0])*10, 10), stoichiometry[0] - stoichiometry[1], stoichiometry[0] + stoichiometry[1], alpha=0.2)

        #plot stoichiometry of different experiments
        for e,experiment in enumerate(intensity_experiments):
            axis[1][0].hist(np.array(experiment[-1], dtype = np.float)/self.calibrations[e]*32, alpha=0.2)
            stoichiometry = np.array(list(zip(*[(np.mean(i),np.std(i)) for i in experiment])))/self.calibrations[0]*32
            axis[1][1].plot(np.arange(0, len(stoichiometry[0])*10, 10), stoichiometry[0], label = experiment_names[e])
            axis[1][1].fill_between(np.arange(0, len(stoichiometry[0])*10, 10), stoichiometry[0] - stoichiometry[1], stoichiometry[0] + stoichiometry[1], alpha=0.2)
            print(e)
        axis[1][1].legend()
        fig.show()

    def analyzeResults(self, calibrated):
        self.update_status_bar("Showing results ...")
        fig, axis = plt.subplots(1,2)
        axis[0].set_ylabel("counts")
        axis[0].set_xlabel("brightness")
        #get data from string to float

        str = self.stch_analysis.export_brightness_3()
        splitstr = str.split("\n")[3].split(",")[0:-2]

        aio = np.array(splitstr).astype(np.float) #all in one
        int_calibration = []
        int_data = []
        for i,replicate in enumerate(self.stch_analysis.replicates):
            single_data = []
            for j in range(len(replicate.particles)-1):
                if replicate.summary['calibration']:
                    int_calibration.append(float(splitstr[0]))
                    splitstr.pop(0)
                else:
                    single_data.append(float(splitstr[0]))
                    splitstr.pop(0)
            if not replicate.summary['calibration']:
                int_data.append(single_data)
        if not calibrated:
            self.stoich_calibration = np.average(np.array(int_calibration).flatten())

        stoichiometry = []
        n_particles = []
        for dat in int_data:
            stoichiometry.append(np.average(np.array(dat))/self.stoich_calibration*32)
            n_particles.append(len(dat))
            axis[0].hist(np.array(dat)/self.stoich_calibration*32, alpha = 0.6)
        axis[1].set_ylabel("stoichiometry")
        axis[1].set_xlabel("time after treatment (min)")
        axis[1].plot(np.arange(0, len(stoichiometry)*10, 10), stoichiometry)
        fig.show()
        return (int_data)

    def onUpdateBrightness(self, event):
        if 0.0 <= self.widgets[self.ba_efficiency.GetId()].getValue() <= 1.0:
            self.update_brightness_tab()
        else:
            wx.MessageBox(
                'Label efficiency should be between 0.0 and 1.0', 'Info', wx.OK | wx.ICON_ERROR)

    def onBrightnessAnalysis(self, event):
        replicate = self.stch_analysis.get_replicate_by_name(
            event.EventObject.Label, load_video=False)
        replicate.summary['ba_analysis'] = event.EventObject.Value

    def update_brightness_tab(self):
        self.update_file_list(self.page_brightness,
                              self.ba_list_spanel,
                              self.ba_list_spanel_sizer,
                              'ba_analysis',
                              self.onBrightnessAnalysis,
                              'no_calibration')
        self.update_brightness_image()

    def update_brightness_image(self):
        self.ba_right_sizer.Hide(self.brightness_panel)
        # self.ba_right_sizer.RemoveChild(self.brightness_panel)
        self.brightness_panel.Destroy()
        self.brightness_panel = wx.Panel(self.page_brightness)
        panel_sizer = wx.BoxSizer(wx.VERTICAL)

        brightness_method = self.ba_brightness_method.GetValue()
        background_method = self.ba_bckg_method.GetValue()
        fitting_method = self.ba_fitting_method.GetValue()
        label_efficiency = float(self.ba_efficiency.GetValue())

        if self.stch_analysis is not None and \
                self.stch_analysis.brightness_analysis(self.methods_val[brightness_method],
                                                       self.methods_val[background_method],
                                                       self.methods_val[fitting_method],
                                                       label_efficiency):
            self.ba_monomers_count.SetLabel(
                str(len(self.stch_analysis.ba_monomers_intensities)))
            self.ba_particles_count.SetLabel(
                str(len(self.stch_analysis.ba_all_intensities)))
            self.ba_monomer_intensity.SetLabel(
                "{:.2f}".format(self.stch_analysis.ba_monomer_intensity))
            self.ba_accuracy_avg.SetLabel(
                "{:.2f}".format(self.stch_analysis.accuracy_avg))

            figure = Figure()
            #axes = figure.add_axes([0.08, 0.08, 0.88, 0.88])
            axes = figure.add_subplot(121)
            axes.set_title('Brightness Analysis')
            axes.set_xlabel('Photon count')
            axes.set_ylabel('Frequency (normalized)')
            # Todo: proper values of tics when numbers are too big
            canvas = FigureCanvas(self.brightness_panel, 1, figure)
            if self.methods_val[fitting_method] == 'gaussfit':
                bar_width = self.stch_analysis.ba_all_intensities_xs[1] - \
                    self.stch_analysis.ba_all_intensities_xs[0]
                axes.bar(self.stch_analysis.ba_all_intensities_xs,
                         self.stch_analysis.ba_all_intensities_ys,
                         bar_width, color='b', alpha=0.1)
            else:
                axes.plot(self.stch_analysis.ba_all_intensities_xs,
                          self.stch_analysis.ba_all_intensities_ys,
                          color='b', alpha=0.3)

            axes.plot(self.stch_analysis.ba_all_intensities_xs,
                      self.stch_analysis.multiple_fitting_res['all'],
                      color='r', alpha=0.7)
            for i, g in enumerate(self.stch_analysis.multiple_fitting_res['ind']):
                axes.plot(self.stch_analysis.ba_all_intensities_xs, g['values'],
                          ls='--', color='r')

            axes = figure.add_subplot(222)
            axes.set_title('Calibration')
            axes.set_xlabel('Photon count')
            axes.set_ylabel('Frequency (normalized)')
            mono_model = self.stch_analysis.mono_model
            mono_xs = np.linspace(
                mono_model['interval'][0], mono_model['interval'][1], 50)
            if self.methods_val[fitting_method] == 'gaussfit':
                bar_width = self.stch_analysis.mono_model['hist_bin_centers'][1] - \
                    self.stch_analysis.mono_model['hist_bin_centers'][0]
                axes.bar(mono_model['hist_bin_centers'],
                         mono_model['hist_values'],
                         bar_width, alpha=0.3, color='b')
                axes.locator_params(nbins=4)
                mono_ys = gaussian(
                    mono_xs, mono_model['mean'], mono_model['sdev'], mono_model['gauss_a'])
            else:
                mono_ys = mono_model['pdf_function'](mono_xs)
            axes.plot(mono_xs, mono_ys, color='r', alpha=0.5)
            axes.axvline(self.stch_analysis.ba_monomer_intensity,
                         ls=':', color='r')

            axes = figure.add_subplot(224)
            axes.set_title('Distribution')
            axes.set_xlabel('Species')
            axes.set_ylabel('Percentage')
            bar_width = 0.8
            xs = np.arange(1, len(self.stch_analysis.ba_distribution)+1)
            axes.set_xticks(xs)
            xs = np.arange(1, len(self.stch_analysis.ba_distribution)+1)-0.4
            axes.bar(xs,
                     self.stch_analysis.ba_distribution, alpha=0.5,
                     yerr=self.stch_analysis.ba_distribution_errs)

            figure.subplots_adjust(left=0.08, right=0.95,
                                   top=0.95, bottom=0.1, hspace=0.3)

            toolbar = NavigationToolbar2Wx(canvas)
            toolbar.Realize()

            panel_sizer.Add(canvas, 1, wx.EXPAND)
            panel_sizer.Add(toolbar, 0)
        else:
            self.ba_monomers_count.SetLabel("-")
            self.ba_particles_count.SetLabel("-")
            self.ba_monomer_intensity.SetLabel("-")

        self.brightness_panel.SetSizer(panel_sizer)
        self.ba_right_sizer.Add(self.brightness_panel, 1, wx.EXPAND)
        self.page_brightness.Layout()

# Tab Stats methods

    def update_stats_tab(self):
        """
        Display summary and stats by file
        :return:
        """
        # Emptying the sizer
        for cb in self.stats_panel_sizer.GetChildren():
            self.stats_panel_sizer.Hide(cb.Window)
            #self.stats_panel.RemoveChild(cb.Window)

        if self.stch_analysis:
            stats_grid = wx.grid.Grid(self.page_stats)
            self.stats_panel_sizer.Add(stats_grid, 1, wx.EXPAND)
            rows = len(self.stch_analysis.replicates)
            stats_grid.CreateGrid(rows, 11)
            stats_grid.SetColLabelValue(0, "File")
            stats_grid.SetColLabelValue(1, "Calib.")
            stats_grid.SetColLabelValue(2, "Parts")
            stats_grid.SetColLabelValue(3, "Mono.")
            stats_grid.SetColLabelValue(4, "Higher")
            stats_grid.SetColLabelValue(5, "Global\nbckg")
            stats_grid.SetColLabelValue(6, "Local\nbckg")
            stats_grid.SetColLabelValue(7, "Global\naccuracy")
            stats_grid.SetColLabelValue(8, "Local\naccuracy")
            stats_grid.SetColLabelValue(9, "Width\navg.")
            stats_grid.SetColLabelValue(10, "Width\nstd.")

            stats_grid.EnableEditing(False)

            # Filling the grid
            for i in range(rows):
                values = {}
                replicate = self.stch_analysis.replicates[i]
                values[0] = basename(replicate.video_name)
                values[1] = replicate.summary['calibration']
                values[2] = len(replicate.particles)
                values[5] = "{:.2f}".format(
                    replicate.summary['global_bckg_avg'])
                values[6] = "{:.2f}".format(
                    replicate.summary['local_bckg_avg'])
                values[7] = "{:.2f}".format(
                    replicate.summary['accuracy_global'])
                values[8] = "{:.2f}".format(
                    replicate.summary['accuracy_local'])
                values[9] = "{:.2f}".format(
                    replicate.summary['fitted_width_avg'])
                values[10] = "{:.2f}".format(
                    replicate.summary['fitted_width_std'])

                values[3] = 0
                values[4] = 0
                ba_method = self.methods_val[self.ba_brightness_method.GetValue(
                )]
                for p in replicate.particles:
                    if p.signals[ba_method].steps_count == 1:  # Is only None when is created
                        values[3] += 1
                    if p.signals[ba_method].steps_count > 1:
                        values[4] += 1

                for k, v in values.items():
                    stats_grid.SetCellValue(i, k, str(v))
                    if values[1]:
                        stats_grid.SetCellBackgroundColour(
                            i, k, wx.Colour(245, 245, 255))
                    if values[2] == 0:
                        stats_grid.SetCellBackgroundColour(
                            i, k, wx.Colour(255, 245, 245))

            stats_grid.AutoSizeColumns()

        self.page_stats.Fit()