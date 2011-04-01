
class ColorPoint:
    """
    Simple color-storage class; stores way-points on a color ramp
    """
    def __init__(self,idx,col,colType):
        # index, X-coordinate, on a palette
        self.idx = idx
        # color; usually an RGBA quad
        self.color = col
        # One of ColorTypes members
        self.colorType = colType
    def __str__(self):
        return "(Index=%d; Color=(%0.3f,%0.3f,%0.3f,%0.3f); ColorType=%d)" % (self.idx, self.color[0], self.color[1], self.color[2], self.color[3], self.colorType)

class ColorTypes:
    """
    Simple enumerated type for internal color formats
    """
    RGBAi = 0
    RGBAf = 1
    HEX6  = 2
    HEX8  = 3

class ColorRamp:
    """
    Model for a simple color ramp

    Typical usage:
    c = ColorRamp(nColors=256)
    
    # must add boundaries
    c.addPoint(0, (0,0,0,0))
    c.addPoint(255, (1,1,1,1))
    
    # now add other colors
    c.addPoint(32, (128,0,0,0), colType=ColorTypes.RGBAi, colScale=255)
    c.addPoint(64, (128,255,255,255), colType=ColorTypes.RGBAi, colScale=255)

    print c

    RGBAf_ColorRamp = c.getRamp()
    RGBAi_ColorRamp = c.getRamp(colType=ColorTypes.RGBAi, colScale=255)
    RGBAhex_ColorRamp = c.getRamp(colType=ColorTypes.HEX6)
    RGBAhex_ColorRamp = c.getRamp(colType=ColorTypes.HEX8)
    
    """
    # we assume linear ramps for now
    LINEAR = 0
    GAUSSIAN = 1
    EXPONENTIAL = 2

    def __init__(self, nColors, *args, **kwargs):
        # size of this ramp
        self.nColors = nColors
        # the list of RGBA float values
        self.ramp = []
        # ordered array of color indices
        self.keys = {}
        # ready to use; boolean; we need at least two
        # color points to define a ramp
        self.ready = False
        #
        if 'handle' in kwargs:
            self.handle = kwargs['handle']

    def __str__(self):
        s = "Ready to use: " + str(self.ready) + "\n"
        s+= "Keys: " + str(self.keys.keys()) + "\n"
        for k in self.keys:
            s += "Color[%d] = %s\n" % (k,self.keys[k])
        s += "ColorRamp with %d colors follows...\n" % self.nColors
        if self.ready:
            s += str(self.getRamp()) + "\n"
        else:
            s += "[]\n"
        return s

    def addColor(self, idx, col, colType=ColorTypes.RGBAf, colScale=1.0):
        # check user input: color location
        if idx<0 or idx>self.nColors-1:
            print "Error: Invalid color idx given.  Valid values for idx are "
            print "Error: [0,%d), but was given was idx=%d" % (self.nColors-1,idx)
            return
        # check user input: duplicate color
        if idx in self.keys:
            print "Error: Invalid index, a color already exists at this point."
            print "Error: Please remove the color at this index before adding."
            return
        # check user input: color format
        if type(col) != ().__class__ or len(col)!=4:
            print "Error: Colors must be spefied as a RGBA tuple with four values."
            print "Error: %s was given instead." % str(col)
            return
        # check user input: color type format
        if colType not in (ColorTypes.RGBAi, ColorTypes.RGBAf):
            print "Error: Color type specification must be either, "
            print "Error: ColorRamp.RGBAi or ColorRamp.RGBAf"
            return

        userCol = None

        # convert color type if needed
        if colType==ColorTypes.RGBAf:
            userCol = col
        elif colType==ColorTypes.RGBAi:
            userCol = map(lambda c: float(c)/float(colScale), col)
            
        # create a ColorPoint and insert it
        self.keys[idx] = ColorPoint(idx, userCol, colType)

        # is this ramp yet good to use?
        self.updateReady()

        # what else do we need to do to modify the model?
        
    def removeColor(self, idx):
        if idx<0 or idx>self.nColors-1:
            print "Error: Tried to remove a color at a place not on this color ramp."
            print "Error: Valid values for the idx are [0,%d) but idx was %d" % (self.nColors,idx)

        # check user input
        if idx not in self.keys: return
        if idx<0 or idx>self.nColors-1: return

        # remove the point
        del self.keys[idx]

        # is this ramp still good to use?
        self.updateReady()


    def updateReady(self):
        # are we ready to use?
        self.ready = (0 in self.keys and self.nColors-1 in self.keys)


    def updateRamp(self, idx=None):
        
        # if idx is specified then it was either added or removed
        # so adjust the ramp about that point
        
        if not self.ready: 
            # no use in updating a ramp w/o proper colors
            print "Msg: This color ramp is not yet ready to use.  Please add"
            print "Msg: at least two colors at the ramp's extreme points 0 and %d" % (self.nColors-1)
            return
        
        # OPTIMIZATION TODO:
        # if idx!=None and idx no in self.keys, then the point
        # was removed, just update around those pts
        # if idx!=None and does exists in self.keys, then they
        # just added this point, so update the pts around it
        
        self.ramp = []
        
        keyList = self.keys.keys()
        keyList.sort()
        keyList.reverse()
        print "KeyList is: " + str(keyList)

        lowerId = keyList.pop()
        while len(keyList)>0:
            upperId = keyList.pop()
            # number of colors in between
            span = abs(upperId-lowerId)
            # get the actual colors
            lowerCol, upperCol = self.keys[lowerId].color, self.keys[upperId].color

            for x in range(span):
                # linear mixing components
                cUpper = float(x) / float(span)
                cLower = 1.0 - cUpper
                
                self.ramp.append((cLower * lowerCol[0] + cUpper * upperCol[0], 
                                  cLower * lowerCol[1] + cUpper * upperCol[1], 
                                  cLower * lowerCol[2] + cUpper * upperCol[2], 
                                  cLower * lowerCol[3] + cUpper * upperCol[3])) 
            lowerId = upperId
        
        # fix the off-by one error
        self.ramp.append(upperCol)

        assert len(self.ramp)==self.nColors, "ColorRamp Logic Error: This ramp supports %d colors ONLY, but %d were found in the ramp." % (self.nColors, len(self.ramp))

    def getRamp(self, colType=ColorTypes.RGBAf, colScale=1.0):
        # update the ramp and return it
        self.updateRamp()

        if colType==ColorTypes.RGBAf:
            if colScale==1.0:
                return self.ramp
        elif colType==ColorTypes.HEX6:
            colScale = 255
            return map(lambda col: "#%02x%02x%02x" % (colScale*col[0],colScale*col[1],colScale*col[2]), self.ramp)
        elif colType==ColorTypes.HEX8:
            colScale = 255
            return map(lambda col: "#%02x%02x%02x%02x" % (colScale*col[0],colScale*col[1],colScale*col[2],colScale*col[3]), self.ramp)

    def toPhotoImageString(self,nRows=1):
        oneLine = "{" + " ".join(self.getRamp(ColorTypes.HEX6)) + "}"
        if nRows==1:
            return oneLine
        else:
            return " ".join([oneLine]*nRows)
            

    # this belongs in the view
    def toPhotoImage(self,nRows=1,root=None):
        try:
            from Tkinter import PhotoImage
        except ImportError, e:
            print "Error: could not import Tk.  No image."
            print "Error: ", e
            return None
        
        img = PhotoImage(width=self.nColors, height=nRows)
        img.put(self.toPhotoImageString(nRows))

        return img


    def getHandle(self):
        return self.handle
    def setHandle(self,handle):
        self.handle = handle


class ColorRampView:
    """
    This handles the UI for mapping the on-screen
    Tcl/Tk UI components to the model, in a usable fashion.
    """

    def __init__(self):
        pass


class ColorRampController:
    """
    Controller

    Takes Model output and makes a UI component
    that allows interaction; feeds results to 
    the controller

    M = ColorRampModel(256)
    ...
    C = ColorRampController(M)
    img = C.toPhotoImage(nRows=100)
    
    can = Canvas(...)
    can.create_image(img, to=(padX,padY), anchor=NW)
    can.setupHandlers(...)
    """

    def __init__(self):
        pass


if __name__=="__main__":
    c = ColorRamp(256)

    # add some colors
    c.addColor(1,(0,0,0,0))
    print c
    c.addColor(2,(1,1,1,1))
    print c
    c.addColor(250,(0.5, 0.5, 0.5, 0.5))
    print c

    # range checking
    c.addColor(-1, (0,0,0,0))
    print c
    c.addColor(256, (1,2,3,4))
    print c

    # color scaling
    c.addColor(45, (128, 255, 64, 32), colType=ColorTypes.RGBAi, colScale=255)
    print c

    # remove a color
    c.removeColor(2)
    print c
    # range checking
    c.removeColor(-1)
    print c
    # range checking
    c.removeColor(2000)
    print c

    # check ready to use
    c.addColor(0, (0,0,0,0))
    print c
    c.addColor(8, (1.0, 0.4, 0.0, 0.0))
    print c
    c.addColor(255, (1,1,1,1))
    print c

    # check ramp types
    d = ColorRamp(32)
    d.addColor(0, (0,0,0,0), colType=ColorTypes.RGBAi, colScale=255)
    d.addColor(31, (255,255,255,255), colType=ColorTypes.RGBAi, colScale=255)
    d.addColor(15, (1.0, 0.0, 0.0, 1.0))
    
    print "Color Ramp as RGAf"
    print d.getRamp()

    print "Color Ramp as HEX6"
    print d.getRamp(ColorTypes.HEX6)

    print "Color Ramp as HEX8"
    print d.getRamp(ColorTypes.HEX8)
    

    print "Does adding/removing a pt screw up the model?"
    f = ColorRamp(360)
    # end pts
    f.addColor(0, (0,0,0,0))
    f.addColor(359, (1,1,1,1))
    print f

    print "Adding a pt"
    f.addColor(7, (1.0, 0.0, 0.5, 0.25))
    print f

    print "Removing a pt"
    f.removeColor(7)
    print f

    print "Add some more colors"
    f.addColor(90, (1.0, 0.0, 0.0, 1.0))
    f.addColor(270, (0.0, 0.0, 1.0, 1.0))
    f.addColor(180, (0.0, 1.0, 0.0, 1.0))

    print "Checking hex8 vlaues"
    print f.getRamp(ColorTypes.HEX8)


    print "To PhotoImage String: nRows=1"
    print f.toPhotoImageString()

    print "To PhotoImage String: nRows=32"
    print f.toPhotoImageString(16)
    
    try:
        from Tkinter import *
        root = Tk()
        padX, padY = 30, 30
        canvas = Canvas(root,height=10+padY,width=360+padX)

        print "Try to make a color ramp image"
        img = f.toPhotoImage(10)

        canvas.create_image((padX/2, padY/2),image=img,anchor=NW)
        canvas.pack()

        root.mainloop()
    except ImportError, e:
        print "WARNING: Tkinter not installed for this Python version."
        print "WARNING: Skipping the Tkinter test"
    
