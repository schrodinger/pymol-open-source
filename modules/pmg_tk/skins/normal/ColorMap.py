
from Tkinter import *

class ColorMap:
    """
    A simple converter from a ColorRamp (a gradient of nColors) to
    a ColorMap (a graph of data + visual way pts)
    """
    def __init__(self,parent,ramp,width=None,height=None,data=None):
        """
        just produces the canvas/image for
        the ColorRamp
        """
        self.parent = parent
        self.ramp = ramp
        self.width = width
        self.height = height
        self.data = data
        # constants
        self.CIRCLE_RADIUS=3
        # container
        self.canvas = None
        self.canvas_ids = []

    def toCanvas(self,parent=None,canvas=None,ramp=None,width=None,height=None):
        """
        get a correct canvas for this color ramp
        """
        print " +ColorMap::toCanvas"
        if None==parent:
            parent=self.parent
        if None==width:
            width=self.width
        if None==height:
            height=self.height
        if None==ramp:
            ramp = self.ramp
        if None==canvas:
            canvas = self.canvas

        r = self.CIRCLE_RADIUS

        ## if self.canvas!=None:
        ##     # empty the canvas
        ##     for x in self.canvas_ids:
        ##         self.canvas.delete(x)
        ##     self.canvas_ids = []
        ## else:
        ##     # create empty canvas
        ##     self.canvas = Canvas(parent,width=width,height=height,bg="white")
        
        # draw colored circles where the users clicked
        for key in ramp.keys.keys():
            curColorPt = ramp.keys[key]

            # get original values
            origX, origY, origColor = curColorPt.idx, curColorPt.color[3], curColorPt.color
        
            # convert to scaled and Tkinter-friendly formats
            # scale X
            x = int(float(origX)/float(ramp.nColors) * width)
            # scale Y
            y = int(height * (1.0-float(origY)))
            # convert the color from RGBA --> HEX6
            colHEX6 = "#%02x%02x%02x" % (origColor[0]*255., origColor[1]*255., origColor[2]*255.)
            # plot the pt
            unique_id = canvas.create_oval((x-r,y-r,x+r,y+r),fill=colHEX6,tags="colorPt")
            self.canvas_ids.append(unique_id)
        print " -ColorMap::toCanvas"            
        return canvas

    def clearCanvas(self):
        for x in self.canvas_ids:
            self.canvas.delete(x)


if __name__=="__main__":

    # test with a ColorRamp
    try:
        root = Tk()
        from ColorRamp import ColorRamp
        c = ColorRamp(256)
        c.addColor(0, (0,0,0,0.7))
        c.addColor(128, (1,0,0,0.5))
        c.addColor(255, (1,1,1,0.2))
        
        m = ColorMap(root, ramp=c, width=256, height=100)
        
        can = m.toCanvas(root)
        can.pack()
        root.mainloop()
    except ImportError:
        print "Warning: Cannot test with a color ramp"

        
