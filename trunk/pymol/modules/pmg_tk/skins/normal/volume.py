
from Tkinter import *
import pymol
from pymol import cmd
import os
import colorsys

#ugh
import time

from hist import VolumeHist
from pymol.colorramping import ColorRamp

import ModalWindow
##############################################################################
# HSV Color Chooser Modal Dialog
#
# this has module scope
selectedColor = None
# this also has module scope
def setColor(event):
    global selectedColor
    w = event.widget
    X, Y = w.canvasx(event.x), w.canvasy(event.y)
    W,H = w.winfo_reqwidth(), w.winfo_reqheight()
    H2 = H/2
    if Y < H2:
        selectedColor = colorsys.hsv_to_rgb(X/W,Y/H2,1.0)
    else:
        selectedColor = colorsys.hsv_to_rgb(X/W,1.0,1.0-(Y-H2)/(H2))

class SimpleColorChooser(ModalWindow.ModalWindow):
    """
    Simple HSV-based ColorChooser Modal Dialog
    """
    def pack(self,x,y):
        padX = padY = 30
        self.frame = Frame(self,width=360+padX,height=300+padY)
        self.canvas = Canvas(self.frame,width=360+padX,height=200+padY)
        self.canvas.pack(side=TOP)
        self.imgName = ""
        if "PYMOL_DATA" in os.environ:
            self.imgName = os.environ["PYMOL_DATA"] + "/pymol/hsv.ppm"
        else:
            self.imgName = "hsv.ppm"
        self.img = PhotoImage(file=self.imgName)
        self.canvas.create_image((padX/2,padY/2),image=self.img,anchor=NW)
        self.canvas.bind("<Button-1>", self.accept)

        self.xText = StringVar()
        self.yText = StringVar()

        # for exact values
        self.xLabel = Label(self.frame, text=self.data_label)
        self.xLabel.pack(side=LEFT)
        self.xEntry = Entry(self.frame,width=5,textvariable=self.xText)#,validate="key",validatecommand=self.xEntryChanged)
        self.xEntry.insert(0, str(self.data_value))
        self.xEntry.pack(side=LEFT)

        self.yLabel = Label(self.frame,text="Opacity [0..1]:")
        self.yLabel.pack(side=LEFT)
        self.yEntry = Entry(self.frame,width=5,textvariable=self.yText)#,validate="key",validatecommand=self.yEntryChanged)
        self.yEntry.insert(0, str(self.alpha_value))
        self.yEntry.pack(side=LEFT)

        # TODO: get this to open right underneath the user's mouse
        self.geometry("+%d+%d" % (x,y))
        self.frame.pack()
    """
    def xEntryChanged(self, *args):
        print "args = ", args
        print "text = ", self.xText.get()
        try:
            value = float(self.xEntry.get())
        except:
            value = self.data_value    
        self.data_value = value        
        return True

    def yEntryChanged(self, *args):
        try:
            value = float(self.yEntry.get())
        except:
            value = self.alpha_value
        self.alpha_value = value
        return True
    """

    def setDataLabel(self, label):
        self.data_label = label

    def setInitialDataAlpha(self, data_value, alpha_value):
        self.data_value = data_value
        self.alpha_value = alpha_value

##############################################################################
# Transfer Function UI
#
class Volume(Frame):
    def __init__(self,app,parent,*kw,**kwargs):
        """
        Volume Transfer Function Control
        """
        # is-a frame; init the base class
        Frame.__init__(self,parent,*kw,**kwargs)

        self.app = app
        self.parent = parent
        self.deferred = 1
        self.shown = False
        self.active = None
        self.active_ramp, self.active_ramp_canvas = None, None
        self.active_map,  self.active_map_canvas  = None, None

        # the three-basic components that make up this UI
        self.object_list = {}
        self.list_box = None

        # constants
        self.padX = 30
        self.padY = 15
        self.COLOR_RAMP_HEIGHT = 8
        self.COLOR_RAMP_WIDTH  = 360+2*self.padX
        self.COLOR_MAP_HEIGHT  = 124
        self.COLOR_MAP_WIDTH   = 360+2*self.padX

        # for events
        self.start_pt = None
        self.user_is_dragging = False

        # ramp list for postponed update
        self.ramp_update = {}

        if 'cmd' in kwargs:
            self.cmd = kwargs["cmd"]
        else:
            self.cmd = cmd

    def deferred_activate(self):
        """
        don't incur the cost of building this frame until someone uses it
        """
        if self.deferred:
            self.deferred = 0
        self.constructMain()

    def constructMain(self):
        """
        construct the main frame
        """
        # update the list with the known volume objects
        self.update_object_list()
        self.update_listbox()
        
        # init the main components from the current selection
        self.update_transferframe()

        # map the events
        self.update_events()
        self.shown = True

    def update_transferframe(self,fast=False):
        self.update_color_map()
        if not fast:
            self.update_color_ramp()
            self.update_events()
            # update the colors
            if self.active!=None:
                if self.active in self.cmd.get_names("public_objects"):
                    self.cmd.volume_color(self.active, self.active_ramp.getRamp())
                    self.cmd.set_volume_ramp(self.active, self.active_ramp.getRampList())
                    self.cmd.recolor()

    def update_is_needed(self):
        # check removed objects
        pubObj = self.cmd.get_names("public_objects")
        for x in self.object_list:
            if x not in pubObj:
                return True
        # check new objects
        for x in self.cmd.get_names("public_objects"):
            if "object:volume"==cmd.get_type(x):
                if x not in self.object_list.keys():
                    return True
        return False

    def update_object_list(self):
        """
        update the internal state list
        """
        # purge list of removed objects
        pubObj = self.cmd.get_names("public_objects")
        object_list_copy = self.object_list.copy()
        for x in object_list_copy:
            if x not in pubObj:
                del self.object_list[x]
        # for all VOLUME type objects not known to the list
        for x in self.cmd.get_names("public_objects"):
            if "object:volume"==cmd.get_type(x):
                if x not in self.object_list.keys():                   
                    if self.cmd.get_volume_is_updated(x) == 0:
                        continue 
                    # make a default pair for this volume object
                    tmpMap = VolumeHist(self.cmd.get_volume_histogram(x,self.cmd),nBins=64)
                    tmpRamp = ColorRamp(360,name=x)
                    self.object_list[x] = (tmpRamp,tmpMap)
                    if self.ramp_update.has_key(x) and self.ramp_update[x]:
                        tmpRamp.addColor(0, (0,0,1,0))
                        tmpRamp.addColor(359, (0,0,1,0))
                        for data, alpha, col, kind in self.ramp_update[x]:
                            self.addWithoutGUINow(x, data, alpha, col, kind=kind)
                        self.ramp_update[x] = []
                        tmpRamp.updateRamp()
                        self.cmd.volume_color(x, tmpRamp.getRamp())
                        self.cmd.set_volume_ramp(x, tmpRamp.getRampList())
                    else:
                        ramp_list = self.cmd.get_volume_ramp(x, self.cmd)
                        if ramp_list:
                            while ramp_list:
                                tmpRamp.addColor(ramp_list[0], 
                                    (ramp_list[1], ramp_list[2], ramp_list[3], ramp_list[4]))
                                ramp_list = ramp_list[5:]
                        else:                    
                            tmpRamp.addColor(0, (0,0,1,0))
                            tmpRamp.addColor(359, (0,0,1,0))
                            tmpRamp.addColor(200, (0.0, 0.0, 1.0, 0.0))
                            tmpRamp.addColor(210, (0.0, 0.8, 1.0, 0.2))
                            tmpRamp.addColor(220, (0.0, 0.0, 1.0, 0.0))
                else:
                    # need to regenerate the histogram
                    (tmpRamp,tmpMap) = self.object_list[x]
                    tmpMap = VolumeHist(self.cmd.get_volume_histogram(x,self.cmd),nBins=64)
                    self.object_list[x] = (tmpRamp,tmpMap)                    
        if len(self.object_list.keys())!=0:
            # guaranteed to exist
            k = self.object_list.keys()[0]
            self.active_ramp, self.active_map = self.object_list[k]
            self.active = k
            self.update_transferframe()

    def update_listbox(self):
        if self.shown:
            self.list_box.delete(0,END)
        else:
            # create the listbox and append the names
            self.list_box = Listbox(self)
        for x in self.object_list.keys():
            self.list_box.insert(END,x)
        self.list_box.grid(row=0, column=2, rowspan=3, sticky=E)
        self.list_box.selection_set(0)

    def update_color_map(self):
        # update the histrogram and points on it
        h,w = self.COLOR_MAP_HEIGHT, self.COLOR_MAP_WIDTH

        # setup the histrogram
        if self.active_map!=None:
            self.active_map_canvas = self.active_map.toCanvas(self,H=h,W=w,padX=self.padX,padY=self.padY)
        else:
            self.active_map_canvas = Canvas(self,height=h,width=w,bg="white")

        # clear content if redrawing
        self.active_map_canvas.delete("colorPt")
        # add color points to the canvas
        if self.active_ramp!=None:
            self.active_ramp.toCanvas(self.active_map_canvas,w,h,self.padX,self.padY)
        self.active_map_canvas.grid(row=0, column=0, columnspan=2, rowspan=2, sticky=NW)

    def update_color_ramp(self):
        h,w = self.COLOR_RAMP_HEIGHT, self.COLOR_RAMP_WIDTH

        # create the color ramp's canvas
        self.color_ramp_canvas = Canvas(self, height=h, width=w)
        # if active ramp, show it; otherwise show blank canvas
        if self.active_ramp!=None:
            self.rampPhotoImg = self.active_ramp.toPhotoImage(nRows=h)
            self.rampImg = self.color_ramp_canvas.create_image((self.padX,0),image=self.rampPhotoImg,anchor=NW)
        # position the canvas on the frame's grid
        self.color_ramp_canvas.grid(row=2, column=0, columnspan=2, sticky=NW)

    def update_events(self):
        #
        # Add a Point -- left mouse CLICK
        # Remove a Point -- middle mouse CLICK or
        #                -- SHIFT-left mouse CLICK
        # Add an iso-surface -- CONTROL-left mouse CLICK
        # Move a Point -- Left mouse DRAG
        # Edit a Point -- Right mouse CLICK on the point
        self.active_map_canvas.bind("<Button-1>", self.startPoint, add=False)
        self.active_map_canvas.bind("<Control-Button-1>", self.addTriplet, add=False)
        self.active_map_canvas.bind("<Shift-Button-1>", self.removePoint, add=False)
        self.active_map_canvas.bind("<Button-2>", self.removePoint, add=False)
        self.active_map_canvas.bind("<ButtonRelease-1>", self.stopPoint)
        self.active_map_canvas.bind("<Button-3>", self.editPoint, add=False)
        self.active_map_canvas.bind("<B1-Motion>", self.mouseMotion)
        if self.list_box:
            self.list_box.bind("<<ListboxSelect>>", self.volumeSelected, add=False)

    def volumeSelected(self,event):
        self.setActive()
        self.update_transferframe()
        
    def setActive(self,event=None):
        """
        sets the following in all at once:
         * active
         * active_ramp
         * active_map
        """
        idx = self.list_box.curselection()
        if idx!=None and len(idx)!=0:
            idx = int(idx[0])
        else:
            return
        self.list_box.selection_set(idx)
        self.active = self.list_box.get(idx)
        self.active_ramp, self.active_map = self.object_list[self.active]

    def mouseMotion(self,event):
        if self.start_pt and self.start_x != None:
            self.user_is_dragging = True
            X = event.widget.canvasx(event.x) - self.padX
            Y = event.widget.canvasy(event.y) - self.padY
            pt = self.start_pt[0]
            new_x = self.active_ramp.checkPoint(pt, self.last_x, X)
            alpha = 1.0 - (float(Y) / float(self.COLOR_MAP_HEIGHT-2*self.padY))
            alpha = self.active_ramp.yToAlpha(alpha)
            if int(self.start_x) == 0:
                new_x = 1
            elif int(self.start_x) == 359:
                new_x = 358
            self.active_ramp.movePoint(pt, new_x, alpha)
            self.update_transferframe(fast=True)
            self.last_x = new_x

    def startPoint(self,event):
        self.start_pt = self.active_map_canvas.find_closest(event.x, event.y) 
        if self.start_pt:
            self.start_x = self.active_ramp.getPoint(self.start_pt[0]) # event.widget.canvasx(event.x) - self.padX
            self.last_x = self.start_x
#        print "MOUSEDOWN(%d,%d)" % (event.widget.canvasx(event.x),
#                                    event.widget.canvasy(event.y))
#        print "ADJUSTED TO (%d,%d)" % (event.widget.canvasx(event.x)-self.padX,
#                                       event.widget.canvasx(event.y)-self.padY)
#
    def stopPoint(self,event):
        
        if not self.user_is_dragging:
            if self.start_pt:
                self.active_ramp.removePoint(self.start_pt[0])
                self.addPoint(event)            
        self.start_pt = None
        self.start_x = None
        self.user_is_dragging = False
        self.update_transferframe()

    def editPoint(self,event):
        # in canvas coords
        X = event.widget.canvasx(event.x)
        Y = event.widget.canvasy(event.y)

        # get the handle
        pt = self.active_map_canvas.find_closest(event.x, event.y)
        if len(pt)==0:
            return
        else:
            pt = pt[0]

        pos = self.active_ramp.getPoint(pt)
        if pos:
            event.x = pos + self.padX

            # update the model
            self.active_ramp.removePoint(pt)
            self.addPoint(event)

            # update the scene
            self.update_transferframe()

    def addWithoutGUI(self,obj,data,alpha,col,kind="triplet"):
        if not self.ramp_update.has_key(obj):
            self.ramp_update[obj] = []
        self.ramp_update[obj].append((data,alpha,col,kind))
 
    def addWithoutGUINow(self,obj,data,alpha,col,kind="triplet"):
        # select the new volume
        #print "--------------------------------------------------------------------------------"
        #print "obj=%s, data=%f, alpha=%d, col=%s, kind=%s" % (obj,data,alpha,str(col),kind)
        #print "--------------------------------------------------------------------------------"
        """
        try:
            if self.list_box==None:
                self.update_object_list()
            self.update_listbox()
            idx = self.list_box.get(0,'end')
            if not len(idx): 
                return
            idx = idx.index(obj)
        except ValueError:
            return

        #print "index in list is: %s for object %d" % (obj,idx)
        # update the list box
        self.list_box.selection_set(idx)
        self.setActive()
        """ 

        ramp, map = self.object_list[obj]

        # plot the data
        X = int((self.COLOR_MAP_WIDTH-2*self.padX) * (data - map.mData) / (map.MData-map.mData))
        if X<0 or X>self.COLOR_MAP_WIDTH-2*self.padX:
            return
        if alpha < 0.0 or alpha > 1.0:
            return
        """
        self.active_ramp.removePoint(200)
        self.active_ramp.removePoint(210)
        self.active_ramp.removePoint(220)
        """
        if kind=="single":
            ramp.addColor(X, (col[0], col[1], col[2], alpha))
        else:
            diff=2
            ramp.addColor(X, (col[0], col[1], col[2], 0))
            ramp.addColor(X+diff, (col[0], col[1], col[2], alpha))
            ramp.addColor(X+2*diff, (col[0], col[1], col[2], 0))

    def addTriplet(self,event):
        event.triplet=True
        self.addPoint(event)

    def addPoint(self,event):
        # ignore clicks when no volume present
        if self.active==None:
            return
      
        global selectedColor

        self.start_pt = None

        # convert the points from CANVAS to COLOR_MAP
        X = event.widget.canvasx(event.x)-self.padX
        Y = event.widget.canvasy(event.y)-self.padY

#        print "Adding point, user COLOR_MAP clicked: %d,%d" % (X,Y)

        # COORDS now relative to COLOR MAP, not CANVAS
        # check bounds
        if X<0 or X>self.COLOR_MAP_WIDTH-2*self.padX:
            return
        if Y<0 or Y>self.COLOR_MAP_HEIGHT-2*self.padY:
            return
        # get color
        chooser = SimpleColorChooser(parent=self.winfo_toplevel(),title="Simple Color Chooser",buttons=False,user_accept=setColor)
        chooser.setDataLabel("Data value [%.3f..%.3f]:" % (self.active_map.mData, self.active_map.MData))

        data = self.active_map.mData + (self.active_map.MData-self.active_map.mData) * \
                                       (float(X) / float(self.COLOR_MAP_WIDTH-2*self.padX))
        alpha = 1.0 - (float(Y) / float(self.COLOR_MAP_HEIGHT-2*self.padY))
        # alpha is now the ratio linearly proportional to its location on the screen
        # now, map to log-scale
        alpha = self.active_ramp.yToAlpha(alpha)
        chooser.setInitialDataAlpha(data, alpha)
        chooser.pack_and_display(event)
        try:
            new_data = float(chooser.xText.get())
        except: 
            new_data = data
        try:
            new_alpha = float(chooser.yText.get())
        except:
            new_alpha = alpha
        X = int((self.COLOR_MAP_WIDTH-2*self.padX) * (new_data - self.active_map.mData) / (self.active_map.MData-self.active_map.mData))
        if X<0 or X>self.COLOR_MAP_WIDTH-2*self.padX:
            return
        alpha = new_alpha        
        if alpha < 0.0 or alpha > 1.0:
            return
        col = selectedColor
        if col!=None:
            if getattr(event,"triplet",None)==None:
                self.active_ramp.addColor(X, (col[0], col[1], col[2], alpha))
            else:
                diff=2
                self.active_ramp.addColor(X, (col[0], col[1], col[2], 0))
                self.active_ramp.addColor(X+diff, (col[0], col[1], col[2], alpha))
                self.active_ramp.addColor(X+2*diff, (col[0], col[1], col[2], 0))
            self.update_transferframe()

    def removePoint(self, event):
        # in canvas coords
        X = event.widget.canvasx(event.x)
        Y = event.widget.canvasy(event.y)

        # get the handle
        pt = self.active_map_canvas.find_closest(X,Y)
        if len(pt)==0:
            return
        else:
            pt = pt[0]

        # update the model
        self.active_ramp.removePoint(pt)

        # update the scene
        self.update_transferframe()


# on the color chooser, add a Entry boxes for density values and opacity values
# test volumes on mobile chipsets and other non-awesome video cards
# save color ramp info to ObjectVolumeState
# 
# (2) edit a pt-color via mouse click
# (7) Users must be able to specify a iso-level
# (5) BUG: load a map X.  Create volume X_volume.  Delete the map.  Load the map.  Color -- CRASH
# (6) save ramps, somehow
# (8) Histogram of data
