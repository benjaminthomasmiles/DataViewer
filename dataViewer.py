
#######################################################
##############VERSION 2.1##############################
##############BMiles-5/3/2013#########################

import numpy as np
from numpy import array, linspace
from scipy import stats
import scipy.io as scio
import gc
import os
import cPickle as pickle
import wx.stc
# Enthought library imports
from traits.api import HasTraits, Instance, Str, Button, Bool, Enum
from traitsui.api import View, Item, UItem, HGroup, VGroup, Group

from pyface.api import FileDialog, OK
from enable.api import ComponentEditor

# Chaco imports
from chaco.api import ArrayPlotData, ColorBar, GridPlotContainer,\
        OverlayPlotContainer, jet, LinearMapper, Plot, GridDataSource,DataRange2D,\
        ImageData, DataRange1D
from chaco.tools.api import LineInspector, PanTool, ZoomTool
from chaco.plot_graphics_context import PlotGraphicsContext

from enthought.chaco.tools.image_inspector_tool import ImageInspectorTool, \
	     ImageInspectorOverlay


class ImagePlot(HasTraits):

    #create UI interface
    plot = Instance(GridPlotContainer)
    LineScans = Instance(GridPlotContainer)
    dataset = {"data1" : "", "data2" : "", "data3" : "", "data4" : "", "notes" : "", "shape" : "", "range" : []}
    value = Str
    saved = Bool(True)
    filedir = Str
    filename = Str
    path = Str
    xLo = Str
    xHi = Str
    yLo = Str
    yHi = Str
    Cmin = Str
    Cmax = Str
    Cmin2 = Str
    Cmax2 = Str
    loaded = False

    Rtrue = Bool(True)
    R2true = Bool(False)
    Fluortrue = Bool(False)
    #Catch variables
    bounds = []
    xLow = []
    yLow = []
    xHigh = []
    yHigh = []

    vscale = None
    notes = Str

    #UI buttons
    reload_button = Button("Reload")
    save_plot_button = Button("Save Plots..")
    save_LineScans_button = Button("Save LineScans..")
    load_button = Button("Load")
    catch_bounds = Button("Catch Bounds")
    to_csv = Button("Generate .CSV")
    scale_set = Button("Set Scale")
    reset = Button("Reset")
    nextfile = Button("Next")
    prevfile = Button("Prev")
    flatten = Button("Flatten")
    readline = Button("LineScan")
    colormap = Button("ApplyColorMap")
    genfilelist = Button("Generate File List")
    out_to_mat = Button("Generate .mat")

    presmooth = Bool
    plot1button = Bool
    plot2button = Bool
    plot3button = Bool
    plot4button = Bool
    fastline = Str
    slowline = Str

    fastscanline = 0
    slowscanline = 0

    fastline  = str(fastscanline)
    slowline  = str(slowscanline)

    fastlinevar = Str
    fastlinemean = Str
    fastlinemax = Str
    fastlinemin = Str
    slowlinevar = Str
    slowlinemean = Str
    slowlinemax = Str
    slowlinemin = Str

    x_f = Bool(label = "Fast Axis x:")
    y_f = Bool(label = "y:")
    z_f = Bool(label = "z:")
    x_s = Bool(label = "Slow Axis x:")
    y_s = Bool(label = "y:")
    z_s = Bool(label = "z:")
    xbounds=None
    ybounds=None

    #wildcard patterns
    file_wildcard = Str("zis File (*.zis)|*.zis|All files|*")
    file_wildcard2 = Str("png File (*.png)|*.png|All files|*")
    file_wildcard3 = Str("mat File (*.mat)|*.mat|All files|*")
    fastlinevars = Group(   Item('fastlinevar'),
                            Item('fastlinemean'),
                            Item('fastlinemax'),
                            Item('fastlinemin'),show_border=True)
    slowlinevars = Group(   Item('slowlinevar'),
                            Item('slowlinemean'),
                            Item('slowlinemax'),
                            Item('slowlinemin'),show_border=True)
    linescangroup = HGroup(Item('LineScans', editor=ComponentEditor(),show_label=False),
                        VGroup(
                            HGroup(
                            Item('plot1button'),
                            Item('plot2button'),
                            ),
                            HGroup(
                            Item('plot3button'),
                            Item('plot4button'),
                            ),
                            HGroup(
                            Item('slowline'),
                            Item('fastline'),),
                            UItem('readline'),
                            fastlinevars,
                            slowlinevars
                            ),label = "LineScan")
    colorgroup = VGroup(HGroup(
                            Item('Cmin',width = -50),
                            Item('Cmax',width = -50),
                            Item('Cmin2',width = -50),
                            Item('Cmax2',width = -50),
                            ),
                        HGroup(
                            UItem('colormap'),
                            UItem('reload_button'),
                            UItem('genfilelist'),
                            UItem('out_to_mat'),
                            ),
                        show_border =True)

    tabs = Group(Item('plot', editor=ComponentEditor(),show_label=False),
                linescangroup,
                Item('notes', style = 'custom', width=.1), layout = "tabbed")
    plotA = Enum("R", "Phase", "R0","R/R0","Fluor", "X", "Y")
    plotB = Enum("R", "Phase", "R0","R/R0","Fluor", "X", "Y")
    plotC = Enum("R", "Phase", "R0","R/R0","Fluor", "X", "Y")
    plotD = Enum("R", "Phase", "R0","R/R0","Fluor", "X", "Y")

    lefttop = Group(
                VGroup(
                    HGroup(
                        Item('value', label = "File", width = 300),
                        Item('presmooth'),
                        ),
                    HGroup(
                        UItem('save_plot_button'),
                        UItem('save_LineScans_button'),
                        UItem('load_button'),
                        UItem('nextfile'),
                        UItem('prevfile'),
                        ),
                    HGroup(
                        Item('Rtrue'),
                        Item('R2true'),
                        Item('Fluortrue'),
                        ),
                    HGroup(
                    Item(name="plotA", style = 'simple'),
                    Item(name="plotB", style = 'simple'),
                    ),
                    HGroup(
                    Item(name="plotC", style = 'simple'),
                    Item(name="plotD", style = 'simple'),
                    ),
                    ),
                )

    righttop = Group(
                    VGroup(
                        HGroup(
                            Item('xLo',width = -50, height = 25),
                            Item('xHi',width = -50),
                            Item('yLo',width = -50, height = 25),
                            Item('yHi',width = -50),
                            ),
                        HGroup(
                            UItem('flatten'),
                            UItem('catch_bounds'),
                            UItem('scale_set'),
                        ),
                    ),
                show_border = True)

    traits_view = View(
            VGroup(
                    HGroup(
                        lefttop,
                        righttop,
                        colorgroup,
                    ),
                    tabs,
            ),
            width=1200, height=700, resizable=True, title="Scan Image Reconstruction")



    def _Rtrue_changed(self):
        if self.Rtrue:
            self.R2true = False
            self.Fluortrue = False
    def _R2true_changed(self):
        if self.R2true:
            self.Rtrue = False
            self.Fluortrue = False
    def _Fluortrue_changed(self):
        if self.Fluortrue:
            self.Rtrue = False
            self.R2true = False
    def _plot1button_changed(self):
        if self.plot1button:
            self.plot2button = False
            self.plot3button = False
            self.plot4button = False
            self._plot2()
    def _plot2button_changed(self):
        if self.plot2button:
            self.plot1button=False
            self.plot3button = False
            self.plot4button = False
            self._plot2()
    def _plot3button_changed(self):
        if self.plot3button:
            self.plot1button=False
            self.plot2button = False
            self.plot4button = False
            self._plot2()
    def _plot4button_changed(self):
        if self.plot4button:
            self.plot1button=False
            self.plot2button = False
            self.plot3button = False
            self._plot2()

    def _colormap_fired(self):
        self.plot1.color_mapper.range.low = float(self.Cmin)
        self.plot1.color_mapper.range.high = float(self.Cmax)

        self.plot2.color_mapper.range.low = float(self.Cmin2)
        self.plot2.color_mapper.range.high = float(self.Cmax2)

        self.plot4.color_mapper.range.low = float(self.Cmin)
        self.plot4.color_mapper.range.high = float(self.Cmax)

        if self.plot1button or self.plot3button or self.plot4button:
            self.plot1lines.color_mapper.range.low = float(self.Cmin)
            self.plot1lines.color_mapper.range.high = float(self.Cmax)
        if self.plot2button:
            self.plot1lines.color_mapper.range.low = float(self.Cmin2)
            self.plot1lines.color_mapper.range.high = float(self.Cmax2)

        #self._plot()
        print "Color Mapped"
    def _scale_set_fired(self):

        self.plot1.range2d.x_range.low = float(self.xLo)
        self.plot1.range2d.x_range.high = float(self.xHi)
        self.plot1.range2d.y_range.low = float(self.yLo)
        self.plot1.range2d.y_range.high = float(self.yHi)

        self.plot2.range2d.x_range.low = float(self.xLo)
        self.plot2.range2d.x_range.high = float(self.xHi)
        self.plot2.range2d.y_range.low = float(self.yLo)
        self.plot2.range2d.y_range.high = float(self.yHi)

        self.plot3.range2d.x_range.low = float(self.xLo)
        self.plot3.range2d.x_range.high = float(self.xHi)
        self.plot3.range2d.y_range.low = float(self.yLo)
        self.plot3.range2d.y_range.high = float(self.yHi)

        self.plot4.range2d.x_range.low = float(self.xLo)
        self.plot4.range2d.x_range.high = float(self.xHi)
        self.plot4.range2d.y_range.low = float(self.yLo)
        self.plot4.range2d.y_range.high = float(self.yHi)

        self.plot1lines.range2d.x_range.low = float(self.xLo)/self.dataset["range"][0][1]*self.dataset["shape"][0]
        self.plot1lines.range2d.x_range.high = float(self.xHi)/self.dataset["range"][0][1]*self.dataset["shape"][0]
        self.plot1lines.range2d.y_range.low = float(self.yLo)/self.dataset["range"][1][1]*self.dataset["shape"][1]
        self.plot1lines.range2d.y_range.high = float(self.yHi)/self.dataset["range"][1][1]*self.dataset["shape"][1]

        self.slow_plot.range2d.x_range.low = float(self.xLo)
        self.slow_plot.range2d.x_range.high = float(self.xHi)
        self.fast_plot.range2d.x_range.low = float(self.yLo)
        self.fast_plot.range2d.x_range.high = float(self.yHi)


    def _reset_fired(self):
        self.vscale = None
        self.Cmin = ""
        self.Cmax = ""
        self._refresh()
    def _value_changed(self):
        #self.saved = False
        self.value = self.value

    def _refresh(self):
        try: self._plot()
        except: print "Option will be applied when plotting"

    def _readline_fired(self):

        #select which line to scan (this should be input as int value of line number until changed)

        self.fastscanline = int(self.fastline)#float(self.slowline)*float(self.dataset["range"][0][1])/self.dataset["shape"][0])
        self.slowscanline = int(self.slowline)#float(self.fastline)*float(self.dataset["range"][1][1])/self.dataset["shape"][1])

        slowstart = int(float(self.xLo)/float(self.dataset["range"][0][1])*self.dataset["shape"][0])
        slowend = int(float(self.xHi)/float(self.dataset["range"][0][1])*self.dataset["shape"][0]-1)
        faststart = int(float(self.yLo)/float(self.dataset["range"][1][1])*self.dataset["shape"][1])
        fastend = int(float(self.yHi)/float(self.dataset["range"][1][1])*(self.dataset["shape"][1])-1)

        fastarray = np.array(self._image_value.data[:,int(self.slowline)])
        slowarray = np.array(self._image_value.data[int(self.fastline),:])

        self.pd.set_data("line_value2",
                                 self._image_value.data[:,self.fastscanline])
        self.pd.set_data("line_value",
                                 self._image_value.data[self.slowscanline,:])

        self.fastlinevar = str(np.std(fastarray))
        self.fastlinemean = str(np.mean(fastarray))
        self.fastlinemax = str(np.amax(fastarray))
        self.fastlinemin = str(np.amin(fastarray))
        self.slowlinevar = str(np.std(slowarray))
        self.slowlinemean = str(np.mean(slowarray))
        self.slowlinemax = str(np.amax(slowarray))
        self.slowlinemin = str(np.amin(slowarray))

        self.slow_plot.title = "Slowline : " + str(self.fastscanline)
        self.fast_plot.title = "Fastline : " + str(self.fastscanline)

    def _flatten_fired(self):
        self._flatten()

    def _out_to_mat_fired(self):
        dialog = FileDialog(default_filename = self.filename+"_MATLAB_Phase"+str(self.dataset["notes"]["start"][2]), action="save as", wildcard=self.file_wildcard3)
        dialog.open()
        if dialog.return_code == OK:
            savefiledir = dialog.directory
            savefilename = dialog.filename
            path = os.path.join(savefiledir, savefilename)
            scio.savemat(path, self.dataset, True)

    def _flatten(self):
        y=0
        x=0
        for line in self.dataset["data1"]:
            x_axis = np.arange(0,int(len(line)))
            y_axis = np.array(line)
            slope,intercept,r_value,p_value,std_err = stats.linregress(x_axis,y_axis)
            x=0
            for point in line:
                self.dataset["data1"][y,x] = point - (slope*x + intercept)
                x+=1
            y+=1

        y=0
        x=0
        for line in self.dataset["data2"]:
            x_axis = np.arange(0,int(len(line)))
            y_axis = np.array(line)
            slope,intercept,r_value,p_value,std_err = stats.linregress(x_axis,y_axis)
            x=0
            for point in line:
                self.dataset["data2"][y,x] = point - (slope*x + intercept)
                x+=1
            y+=1
        y=0
        x=0
        for line in self.dataset["data3"]:
            x_axis = np.arange(0,int(len(line)))
            y_axis = np.array(line)
            slope,intercept,r_value,p_value,std_err = stats.linregress(x_axis,y_axis)
            x=0
            for point in line:
                self.dataset["data3"][y,x] = point - (slope*x + intercept)
                x+=1
            y+=1
        y=0
        x=0
        for line in self.dataset["data4"]:
            x_axis = np.arange(0,int(len(line)))
            y_axis = np.array(line)
            slope,intercept,r_value,p_value,std_err = stats.linregress(x_axis,y_axis)
            x=0
            for point in line:
                self.dataset["data4"][y,x] = point - (slope*x + intercept)
                x+=1
            y+=1
        self._plot()
        print "Flattened"



    def _catch_bounds_fired(self):
        try:
            self.xLow = float(self.plot1.range2d.x_range.low)
            self.xHigh = float(self.plot1.range2d.x_range.high)
            self.yLow = float(self.plot1.range2d.y_range.low)
            self.yHigh = float(self.plot1.range2d.y_range.high)

            self.xLo = str(self.xLow)
            self.xHi = str(self.xHigh)
            self.yLo = str(self.yLow)
            self.yHi = str(self.yHigh)
        except: print "Please plot first"

    def _save_plot_button_fired(self):
        dialog = FileDialog(default_filename = self.filename+"_Plots_", action="save as", wildcard=self.file_wildcard2)
        dialog.open()
        if dialog.return_code == OK:
            savefiledir = dialog.directory
            savefilename = dialog.filename
            path = os.path.join(savefiledir, savefilename)
            #self.plot.do_layout(force=True)
            plot_gc = PlotGraphicsContext(self.plot.outer_bounds)
            plot_gc.render_component(self.plot)
            plot_gc.save(path)

    def _save_LineScans_button_fired(self):
        dialog = FileDialog(default_filename = self.filename+"_LineScan_", action="save as", wildcard=self.file_wildcard2)
        dialog.open()
        if dialog.return_code == OK:
            savefiledir = dialog.directory
            savefilename = dialog.filename
            path = os.path.join(savefiledir, savefilename)
            self.LineScans.do_layout(force=True)
            plot_gc = PlotGraphicsContext(self.LineScans.outer_bounds)
            plot_gc.render_component(self.LineScans)
            plot_gc.save(path)




    def _load_button_fired(self):
        dialog = FileDialog(action="open", wildcard=self.file_wildcard)
        dialog.open()
        if dialog.return_code == OK:
            self.value = dialog.filename
            title = self.value
            self.path = dialog.path
            self.filedir = dialog.directory
            self.allfiles =[]

            for filenames in os.walk(self.filedir):
                for files in filenames:
                    for afile in files:
                        if ".zis" in str(afile):
                            if ".png" not in str(afile)and ".mat" not in str(afile):
                                self.allfiles.append(self.filedir+"\\"+afile)
            self.filename = dialog.filename
            self.saved = True
            self.loaded = True
            self._plotter(self.path)
            self._plot()

    def _nextfile_fired(self):
        if self.loaded:
            nextone = False
            path2=self.path
##            self.allfiles =[]
##            for filenames in os.walk(self.filedir):
##                for files in filenames:
##                    for afile in files:
##                        if ".zis" in str(afile):
##                            if ".png" not in str(afile):
##                                self.allfiles.append(self.filedir+"\\"+afile)
            for afile in self.allfiles:
                if nextone == True:
                    self.path = afile
                    junk,self.value = afile.split(self.filedir+"\\")
                    self.filename =self.value
                    self._plotter(self.path)
                    self._plot()
                    nextone=False
                    break
                if afile == path2:
                    nextone = True


    def _prevfile_fired(self):
        if self.loaded:
            nextone = False
            path2=self.path
            for afile in self.allfiles:
                if afile == path2:
                    self.path = prevfile
                    junk,self.value = prevfile.split(self.filedir+"\\")
                    self.filename = self.value
                    self._plotter(self.path)
                    self._plot()
                    break
                prevfile = afile

    def _genfilelist_fired(self):
        if self.loaded:
            event = {'trial': 0 , "settings" : "", "notes": "", "time": ""}
            eventlist = {"Description": "", "0" : event}
            i=1
            for afile in self.allfiles:
                #grab file name
                junk,currentfilename = afile.split(self.filedir+"\\")
                print "Working on file : " + currentfilename
                #unpickle file and grab data
                try:
                    currentfile = open(afile,'rb')
                    data = pickle.load(currentfile)
                    currentfile.close()

                    foldername,time = currentfilename.split("_")
                    time,junk = time.split(".")
                    settings = data['settings']['scan']
                    strsettings = ""
                    for key, value in settings.iteritems() :
                        strsettings += str(key) + " " + str(value)+ "\n"
                    newtrial = {'trial': i, "settings" : strsettings, "notes": "", "time": time}
                    eventlist[str(i)] = newtrial

                    i +=1
                except:
                    print "\tcorrupt file, skipping"
            settings = ""
            strsettings =""
            newtrial = ""

            #save data to data logger compatible file
            a = os.getcwd() + "/eventlist_"+foldername
            if not os.path.isdir(a):
                print "made"
                os.makedirs(a)
            b = a+ "/filelist.log"
            f1 = open(b, "w")
            pickle.dump(eventlist, f1)
            f1.close()

            print "File Write Complete"
        else:
            print "Please load a folder first"

    def _reload_button_fired(self):
        self._plotter(self.path)
        self._plot()
    #-----------------------------------------------
    # Private API
    #-----------------------------------------------

    def _plot(self):
        self.notes = ""
        for key, value in self.dataset["notes"].iteritems() :
            self.notes += str(key) + " " + str(value)+ "\n"
        self.container1 = GridPlotContainer(shape = (2,4), spacing = (0,0), use_backbuffer=True,
	                                     valign = 'top', bgcolor = 'white')

        self.plotdata1 = ArrayPlotData(imagedata = self.dataset["data1"])
        self.plotdata2 = ArrayPlotData(imagedata = self.dataset["data2"])
        self.plotdata3 = ArrayPlotData(imagedata = self.dataset["data3"])
        self.plotdata4 = ArrayPlotData(imagedata = self.dataset["data4"])

        s = ""
        f = ""
        if self.x_s: s = "X"
        if self.y_s: s = "Y"
        if self.z_s: s = "Z"
        if self.x_f: f = "X"
        if self.y_f: f = "Y"
        if self.z_f: f = "Z"

        if  self.Rtrue:
            self.plot1 = Plot(self.plotdata1, title = "R : "+ str(self.dataset["notes"]["start"]))
            self.plot2 = Plot(self.plotdata2, title = "PHASE")
            self.plot3 = Plot(self.plotdata3, title = "X")
            self.plot4 = Plot(self.plotdata4, title = "Y" + str(self.dataset["notes"]["start"]))
        if self.R2true:
            self.plot1 = Plot(self.plotdata1, title = "R : "+ str(self.dataset["notes"]["start"]))
            self.plot2 = Plot(self.plotdata2, title = "PHASE")
            self.plot3 = Plot(self.plotdata3, title = "R0")
            self.plot4 = Plot(self.plotdata4, title = "R/R0" + str(self.dataset["notes"]["start"]))
        if self.Fluortrue:
            self.plot1 = Plot(self.plotdata1, title = "R : "+ str(self.dataset["notes"]["start"]))
            self.plot2 = Plot(self.plotdata2, title = "PHASE")
            self.plot3 = Plot(self.plotdata3, title = "Fluor")
            self.plot4 = Plot(self.plotdata4, title = "Y" + str(self.dataset["notes"]["start"]))

        self.plot1.img_plot("imagedata", xbounds = (self.dataset["range"][0][0],self.dataset["range"][0][1]), ybounds = (self.dataset["range"][1][0],self.dataset["range"][1][1]))
        self.plot2.img_plot("imagedata", xbounds = (self.dataset["range"][0][0],self.dataset["range"][0][1]), ybounds = (self.dataset["range"][1][0],self.dataset["range"][1][1]))
        self.plot3.img_plot("imagedata", xbounds = (self.dataset["range"][0][0],self.dataset["range"][0][1]), ybounds = (self.dataset["range"][1][0],self.dataset["range"][1][1]))
        self.plot4.img_plot("imagedata", xbounds = (self.dataset["range"][0][0],self.dataset["range"][0][1]), ybounds = (self.dataset["range"][1][0],self.dataset["range"][1][1]))



#        self.scale = Str(self.plot3.color_mapper.high)
#        plot1.index_axis.title = str(f) + ' (um)'

##ADD TOOLS
        self. plot1.value_axis.title = str(s) + ' (um)'
        self.plot1.tools.append(PanTool(self.plot1))
        zoom1 = ZoomTool(component=self.plot1, tool_mode="box", always_on=False)
        self.plot1.overlays.append(zoom1)

        self.plot2.value_axis.title = str(s) + ' (um)'
        self.plot2.tools.append(PanTool(self.plot2))
        zoom2 = ZoomTool(component=self.plot2, tool_mode="box", always_on=False)
        self.plot2.overlays.append(zoom2)

        self.plot3.index_axis.title = str(f) + ' (um)'
        self.plot3.value_axis.title = str(s) + ' (um)'
        self.plot3.tools.append(PanTool(self.plot3))
        zoom3 = ZoomTool(component=self.plot3, tool_mode="box", always_on=False)
        self.plot3.overlays.append(zoom3)

        self.plot4.index_axis.title = str(f) + ' (um)'
        self.plot4.value_axis.title = str(s) + ' (um)'
        self.plot4.tools.append(PanTool(self.plot4))
        zoom4 = ZoomTool(component=self.plot4, tool_mode="box", always_on=False)
        self.plot4.overlays.append(zoom4)

##ADD COLORBARS
        self.colorbar1 = ColorBar(index_mapper=LinearMapper(range=self.plot1.color_mapper.range),
                        color_mapper=self.plot1.color_mapper,
                        orientation='v',
                        resizable='v',
                        width=20,
                        padding=5)
        self.colorbar1.plot = self.plot1
        self.colorbar1.padding_left = 45
        self.colorbar1.padding_right= 5

        self.colorbar2 = ColorBar(index_mapper=LinearMapper(range=self.plot2.color_mapper.range),
                        color_mapper=self.plot2.color_mapper,
                        orientation='v',
                        resizable='v',
                        width=20,
                        padding=5)
        self.colorbar2.plot = self.plot2
        self.colorbar2.padding_left = 10

        self.colorbar3 = ColorBar(index_mapper=LinearMapper(range=self.plot3.color_mapper.range),
                        color_mapper=self.plot3.color_mapper,
                        orientation='v',
                        resizable='v',
                        width=20,
                        padding=5)
        self.colorbar3.plot = self.plot3
        self.colorbar3.padding_left = 45
        self.colorbar3.padding_right= 5

        self.colorbar4 = ColorBar(index_mapper=LinearMapper(range=self.plot4.color_mapper.range),
                        color_mapper=self.plot4.color_mapper,
                        orientation='v',
                        resizable='v',
                        width=20,
                        padding=5)
        self.colorbar4.plot = self.plot4
        self.colorbar4.padding_left = 15



        self.colorbar1.tools.append(PanTool(self.colorbar1, constrain_direction="y", constrain=True))
        self.zoom_overlay1 = ZoomTool(self.colorbar1, axis="index", tool_mode="range",
                            always_on=True, drag_button="right")
        self.colorbar1.overlays.append(self.zoom_overlay1)

        self.colorbar2.tools.append(PanTool(self.colorbar2, constrain_direction="y", constrain=True))
        self.zoom_overlay2 = ZoomTool(self.colorbar2, axis="index", tool_mode="range",
                            always_on=True, drag_button="right")
        self.colorbar2.overlays.append(self.zoom_overlay2)

        self.colorbar3.tools.append(PanTool(self.colorbar3, constrain_direction="y", constrain=True))
        self.zoom_overlay3 = ZoomTool(self.colorbar3, axis="index", tool_mode="range",
                            always_on=True, drag_button="right")
        self.colorbar3.overlays.append(self.zoom_overlay3)

        self.colorbar4.tools.append(PanTool(self.colorbar4, constrain_direction="y", constrain=True))
        self.zoom_overlay4 = ZoomTool(self.colorbar4, axis="index", tool_mode="range",
                            always_on=True, drag_button="right")
        self.colorbar4.overlays.append(self.zoom_overlay4)




        self.container1.add(self.colorbar1)
        self.container1.add(self.plot1)
        self.container1.add(self.plot2)
        self.container1.add(self.colorbar2)
        self.container1.add(self.colorbar3)
        self.container1.add(self.plot3)
        self.container1.add(self.plot4)
        self.container1.add(self.colorbar4)

        self.plot1.padding_right = 5
        self.plot2.padding_left = 5
        self.plot1.padding_bottom = 15
        self.plot2.padding_bottom = 15
        self.plot3.padding_top = 15
        self.plot4.padding_top = 15
        self.plot1.x_axis.orientation = "top"
        self.plot2.x_axis.orientation = "top"
        self.plot2.y_axis.orientation = "right"
        self.plot3.padding_right = 5
        self.plot4.padding_left = 5
        self.plot4.y_axis.orientation = "right"

        self.colorbar1.padding_top = self.plot1.padding_top
        self.colorbar1.padding_bottom = self.plot1.padding_bottom
        self.colorbar2.padding_top = self.plot2.padding_top
        self.colorbar2.padding_bottom = self.plot2.padding_bottom
        self.colorbar3.padding_top = self.plot3.padding_top
        self.colorbar3.padding_bottom = self.plot3.padding_bottom
        self.colorbar4.padding_top = self.plot4.padding_top
        self.colorbar4.padding_bottom = self.plot4.padding_bottom

        imgtool = ImageInspectorTool(self.plot1)
        self.plot1.tools.append(imgtool)
        overlay = ImageInspectorOverlay(component=self.plot1, image_inspector=imgtool,
                                    bgcolor="white", border_visible=True)
        self.plot1.overlays.append(overlay)

        if not(self.plot1button or self.plot2button or self.plot3button or self.plot4button):
            self.plot1button = True

        self.plot = self.container1

        self._plot2()




#######line plots##############################################################################################

    def _plot2(self):

        if self.plot1button:
            self.lineplotdata = self.dataset["data1"]
            title = 'R : ' + str(self.dataset["notes"]["start"])
        if self.plot2button:
            self.lineplotdata = self.dataset["data2"]
            title = 'Phase : ' + str(self.dataset["notes"]["start"])
        if self.plot3button:
            self.lineplotdata = self.dataset["data3"]
            title = 'R0 : ' + str(self.dataset["notes"]["start"])
        if self.plot4button:
            self.lineplotdata = self.dataset["data4"]
            title = 'R/R0 : ' + str(self.dataset["notes"]["start"])

        self.plotdata5 = ArrayPlotData(imagedata = self.lineplotdata)



        self._image_index = GridDataSource(array([]),
                                          array([]),
                                          sort_order=("ascending","ascending"))

        self.xs = linspace(self.dataset["range"][0][0], self.dataset["range"][0][1], self.dataset["shape"][0])
        self.ys = linspace(self.dataset["range"][1][0], self.dataset["range"][1][1], self.dataset["shape"][1])

        self._image_index.set_data(self.xs, self.ys)

        image_index_range = DataRange2D(self._image_index)


        self._image_value = ImageData(data=array([]), value_depth=1)
        self._image_value.data = self.lineplotdata
        image_value_range = DataRange1D(self._image_value)

        s = ""
        f = ""
        if self.x_s: s = "X"
        if self.y_s: s = "Y"
        if self.z_s: s = "Z"
        if self.x_f: f = "X"
        if self.y_f: f = "Y"
        if self.z_f: f = "Z"


        self.plot1lines = Plot(self.plotdata5,  title = title)
        self.plot1lines.img_plot("imagedata", xbounds = (self.dataset["range"][0][0],self.dataset["range"][0][1]), ybounds = (self.dataset["range"][1][0],self.dataset["range"][1][1]), colormap=jet)
        img_plot = self.plot1lines.img_plot("imagedata")[0]
        imgtool = ImageInspectorTool(img_plot)
        img_plot.tools.append(imgtool)
        overlay = ImageInspectorOverlay(component=img_plot, image_inspector=imgtool,
                                    bgcolor="white", border_visible=True)
        self.plot1lines.overlays.append(overlay)
##ADD TOOLS

        self.plot1lines.tools.append(PanTool(self.plot1lines))
        zoom1 = ZoomTool(component=self.plot1lines, tool_mode="box", always_on=False)
        self.plot1lines.overlays.append(zoom1)

        self.plot1lines.overlays.append(LineInspector(component=self.plot1lines,
                                               axis='index_x',
                                               inspect_mode="indexed",
                                               write_metadata=True,
                                               is_listener=False,
                                               #constrain_key="right",
                                               color="white"))
        self.plot1lines.overlays.append(LineInspector(component=self.plot1lines,
                                               axis='index_y',
                                               inspect_mode="indexed",
                                               write_metadata=True,
                                               color="white",
                                               is_listener=False))




##ADD COLORBAR

        self.colorbar5  = ColorBar(index_mapper=LinearMapper(range=self.plot1lines.color_mapper.range),
                        color_mapper=self.plot1lines.color_mapper,
                        orientation='v',
                        resizable='v',
                        width=20,
                        padding=5)
        self.colorbar5.plot = self.plot1lines


        self.colorbar5.tools.append(PanTool(self.colorbar5, constrain_direction="y", constrain=True))
        self.zoom_overlay5 = ZoomTool(self.colorbar5, axis="index", tool_mode="range",
                            always_on=True, drag_button="right")
        self.colorbar5.overlays.append(self.zoom_overlay5)

        self.pd = ArrayPlotData(line_index = array([]),
                                line_value = array([]))



        self.slow_plot = Plot(self.pd,   title = "Slowline : " + self.slowline)
        self.slow_plot.plot(("line_index", "line_value"),
                             line_style='solid')

        self.pd.set_data("line_index2", array([]))
        self.pd.set_data("line_value2", array([]))

        self.fast_plot = Plot(self.pd,   title = "Fastline : " + self.fastline)
        self.fast_plot.plot(("line_index2", "line_value2"),
                             line_style='solid')


        self.pd.set_data("line_index", self.xs)
        self.pd.set_data("line_index2", self.ys)
        self.pd.set_data("line_value",
                                 self._image_value.data[self.fastscanline,:])
        self.pd.set_data("line_value2",
                                 self._image_value.data[:,self.slowscanline])

        self.colorbar5.padding= 0
        self.colorbar5.padding_left = 15
        #self.colorbar5.height = 400
        self.colorbar5.padding_top =50
        self.colorbar5.padding_bottom = 0
        self.colorbar5.padding_right = 25
        self.colorbar5.padding_left = 50


        self.plot1lines.width = 300
        self.plot1lines.padding_top = 50

        self.plot1lines.index_axis.title = 'fast axis (um)'
        self.plot1lines.value_axis.title = 'slow axis (um)'


        self.slow_plot.width = 100
        self.slow_plot.padding_right = 20


        self.fast_plot.width = 100
        self.fast_plot.padding_right = 20



        self.container2 = GridPlotContainer(shape = (1,2), spacing = ((0,0)), use_backbuffer=True,
	                                     valign = 'top', halign = 'center', bgcolor = 'white')
        self.container3 = GridPlotContainer(shape = (2,1), spacing = (0,0), use_backbuffer=True,
	                                     valign = 'top', halign = 'center', bgcolor = 'grey')
        self.container4 = GridPlotContainer(shape = (1,2), spacing = (0,0), use_backbuffer=True,
	                                     valign = 'top', halign = 'center', bgcolor = 'grey')

        self.container2.add(self.colorbar5)
        self.container3.add(self.fast_plot)
        self.container3.add(self.slow_plot)
        self.container4.add(self.container3)
        self.container4.add(self.plot1lines)
        self.container2.add(self.container4)

        self.LineScans = self.container2

        self._readline_fired()
        self._scale_set_fired()

    def _load_file(self,i):
        file = open(str(i),'rb')
        print "\n\n" + str(i)   + " is open"
        data = pickle.load(file)
        file.close()
        print "\nfile.closed"
        return data

    def _pick_plots(self):
        #plotA = Enum("R", "Phase", "R0","R/R0","Fluor", "X", "Y")

        if self.dataset['plotA']== "R":
            print "hello you called?"

    def _plotter(self, i):
        #open file with specified path unpickle and prepare to pass to dictionary
        file = open(str(i),'rb')
        print "\n\n" + str(i)   + " is open"
        self.data = pickle.load(file)

##        for key, value in self.data['data'].iteritems():
##            print key
        file.close()
        print "\nfile.closed"



        #--APPLY LOGIC CONDITIONS TO DATA AND ASSIGN TO PLOTING 2D ARRAYS

        #load ->

        self.dataset.clear()
        self.dataset = {"plota": self.plotA, "data1" : "", "data2" : "", "data3" : "", "data4" : "", "notes" : self.data['settings']['scan'], "shape" : "", "range" : []}
        self.dataset["range"] = [[0,15],[0,15]]
        gc.collect()

        self._pick_plots()
        print "Display is %d by %d um, and %d by %d by %d datapoints" %(self.dataset["notes"]['range'][0], self.dataset["notes"]['range'][1],self.dataset["notes"]['npoints'][0],self.dataset["notes"]['npoints'][1],self.dataset["notes"]['npoints'][2])
        if self.dataset["notes"]['axes'][0] == 0:
            fast = self.dataset["notes"]['npoints'][0]
            self.dataset["range"][0][1] = self.dataset["notes"]["fast_axis_range"]
        elif self.dataset["notes"]['axes'][0] == 1:
            slow = self.dataset["notes"]['npoints'][0]
            self.dataset["range"][1][1] = self.dataset["notes"]['range'][0]
        if self.dataset["notes"]['axes'][1] == 0:
            fast = self.dataset["notes"]['npoints'][1]
            self.dataset["range"][0][1] = self.dataset["notes"]["fast_axis_range"]
        elif self.dataset["notes"]['axes'][1] == 1:
            slow = self.dataset["notes"]['npoints'][1]
            self.dataset["range"][1][1] = self.dataset["notes"]['range'][1]
        if self.dataset["notes"]['axes'][2] == 0:
            fast = self.dataset["notes"]['npoints'][2]
            self.dataset["range"][0][1] = self.dataset["notes"]["fast_axis_range"]
        elif self.dataset["notes"]['axes'][2] == 1:
            slow = self.dataset["notes"]['npoints'][2]
            self.dataset["range"][1][1] = self.dataset["notes"]['range'][2]

        if self.xLo == "":
            self.xLo = str(self.dataset["range"][0][0])
            self.xHi = str(self.dataset["range"][0][1])
            self.yLo = str(self.dataset["range"][1][0])
            self.yHi = str(self.dataset["range"][1][1])

        self.dataset["shape"] =(fast,slow)
        dataX   = np.array(self.data['data']['lia_x']["0"])
        dataY   = np.array(self.data['data']['lia_y']["0"])
        if self.Fluortrue:
            data3   = np.array(self.data['data']['dio2'])
        else:
            data3 = np.array(self.data['data']['lia_y']["0"])
        data4   = np.array(self.data['data']['lia_x']["0"])#['auxin0'])


        dataX.shape = (self.dataset["shape"])
        dataY.shape = (self.dataset["shape"])
        data3.shape = (self.dataset["shape"])
        data4.shape = (self.dataset["shape"])

        #flatten X and Y
        if self.presmooth:
            y=0
            x=0
            for line in dataX:
                x_axis = np.arange(0,int(len(line)))
                y_axis = np.array(line)
                slope,intercept,r_value,p_value,std_err = stats.linregress(x_axis,y_axis)
                x=0
                for point in line:
                    dataX[y,x] = point - (slope*x + intercept)
                    x+=1
                y+=1
            y=0
            x=0
            for line in dataY:
                x_axis = np.arange(0,int(len(line)))
                y_axis = np.array(line)
                slope,intercept,r_value,p_value,std_err = stats.linregress(x_axis,y_axis)
                x=0
                for point in line:
                    dataY[y,x] = point - (slope*x + intercept)
                    x+=1
                y+=1

        data1 = np.sqrt(np.multiply(dataX,dataX)+np.multiply(dataY,dataY))
        data2 = np.arctan2(dataY,dataX)
        if self.R2true:
            try:

                dataX2   = np.array(self.data['data']['lia_x']["3"])
                dataY2   = np.array(self.ata['data']['lia_y']["3"])

                dataX2.shape = (self.dataset["shape"])
                dataY2.shape = (self.dataset["shape"])
                data2 = data3
                data3 = np.sqrt(np.multiply(dataX2,dataX2)+np.multiply(dataY2,dataY2))
                data4 = data1/data3
                print "plot1: R \tplot2: Phase"
                print "plot3: R0 \toplot3: R/R0"
            except:
                self.Rtrue = True
                self.R2true = False
                print "No second channel"

        self.dataset["data1"] = data1
        self.dataset["data2"] = data2
        self.dataset["data3"] = data3
        self.dataset["data4"] = data4


        data1=[]
        data2=[]
        data3=[]
        data4=[]
        dataX = []
        dataY = []
        self.data = []
        gc.collect()
        return self.dataset
    #axes, bounds,R, PHASE, AUX, FLUO, notes



def mainPlotter():
    ImagePlot().configure_traits()

    return



mainPlotter()