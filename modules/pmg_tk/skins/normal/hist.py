from Tkinter import *
from random import gauss
import math
import sys

class VolumeHist:

    def __init__(self,hist,nBins=64,root=None,W=None,H=None,padX=None,padY=None):
        """
        hist, nBins, minD, maxD, meanD, stdevD required at creation.  Width, Height, padding 
        and root can be specified at image creation time

        see bottom of file for usage/testing

        NOTE: I think the sd-lines are off by a little, but cannot find why; also SD's only
        relevant for normal distributions...

        TODO: Maybe instead of just data, pass in the filename and get the data from the C-layer?

        """
        self.W = W
        self.H = H
        self.padX = padX
        self.padY = padY
        self.root=root
        self.nBins = nBins
        
        # array of N-bins, all counts set to 0
        self.hist = hist[4:]
        self.hist.append(0.0)
        self.IDstore = []
        
        # min, max, and range
        self.mData = hist[0]
        self.MData = hist[1]
        self.meanData = hist[2]
        self.sdData = hist[3]

        #print "PYTHON HISTOGRAM ", self.mData, self.MData, self.meanData, self.sdData

        self.rangeData = self.MData-self.mData
        self.binStep = float(self.rangeData) / float(self.nBins)

        #print "Histogram", len(self.hist)
        #for x in self.hist:
        #    print x

    def rawToX(self,raw):
        """
        given the raw X value, returns the position on the X
        axis, assume 0 padding and unit width
        """
        return float(raw-self.mData) / self.rangeData

    def sdLevelToX(self,level):
        """
        given the SD level, returns the X-position for this point
        on the X axis, assuming 0 padding and unit width
        """
        return self.rawToX(self.sdData * float(level) + self.meanData)

    def binToX(self,bin):
        return float(bin) / float(self.nBins)

    def toCanvas(self,root=None,W=None,H=None,padX=None,padY=None):
        """
        do the plotting; returns a canvas with root specified here or
        at initialization time
        """
        if W==None:
            if self.W==None:
                W = 200
            else:
                W = self.W
        if H==None:
            if self.H==None:
                H = 100
            else:
                H = self.H
        if padX==None:
            if self.padX == None:
                padX = 15
            else:
                padX = self.padX
        if padY==None:
            if self.padY == None:
                padY = 10
            else:
                padY = self.padY
        if root==None:
            if self.root==None:
                root = Tk()
            else:
                root = self.root

        # dimensions of the panel - padding at top and bottom
        w = W-padX*2
        h = H-padY*2

        # Frame/Canvas
        self.c = Canvas(root, height=H, width=W, bg="white")
        self.c.pack(anchor=NW)

        # hist values
        m = min(self.hist)
        M = max(self.hist)
        mean = float(sum(self.hist))/float(len(self.hist))

        # canvas coordos
        left = (padX,)
        right = (w+padX,)
        top = (padY,)
        bottom = (h+padY,)
        topLeft = left+top
        bottomLeft = left+bottom
        topRight = right+top
        bottomRight = right+bottom
        middleX = (padX+(right[0]-left[0])/2,)

        # y-axis
        self.IDstore.append(self.c.create_line(topLeft+bottomLeft, fill="black", width=3,tags="vol_hist"))
        # x-axis
        self.IDstore.append(self.c.create_line(bottomLeft+bottomRight, fill="black", width=3,tags="vol_hist"))
        # axis labels
        R = 10
        fill="black"
        # this skipping stuff is for levels that are too close to print easily; institutes
        # a level padding of ~40 pixels per level
        xSkip = abs(int(40. / (w*(self.sdLevelToX(0) - self.sdLevelToX(1)))))
        if xSkip > 1:
            R *= xSkip
        xSkipCount = 0
        for sd in range(-R,R+1):
            if xSkip>1:
                xSkipCount+=1
                if xSkipCount % xSkip!=0:
                    continue
            xPos = padX + w*self.sdLevelToX(sd)
            if xPos>padX and xPos<w+padX:
                self.IDstore.append(self.c.create_line((xPos,top[0],xPos,bottom[0]), fill="#883333", dash=(2,5),tags="vol_hist"))

                # SIGMA levels
                if sd<0:
                    fill="red"
                else:
                    fill="black"
                self.IDstore.append(self.c.create_text((xPos,padY/2), text=u"%d\u03c3"%sd,tags="vol_hist",font=("helvetica", 9, "normal"), fill=fill))
                # actual values
                if self.meanData+self.sdData*sd<0:
                    fill="red"
                else:
                    fill="black"

                self.IDstore.append(self.c.create_text((xPos,H-padY/2), text=u"%.2f"% (self.meanData+sd*self.sdData),tags="vol_hist",font=("helvetica", 8, "normal"),fill=fill))

        # create a 0 position marker
        xPos = padX + w*self.rawToX(0.0)
        if xPos>padX and xPos<w+padX:
            self.IDstore.append(self.c.create_line((xPos,top[0],xPos,bottom[0]),fill="#ff0000",dash=(2,5),tags="vol_hist"))

        # horizontal tickmarks
        nTicks = 4
        nBase = 10.
        for x in range(nTicks+1):
            yVal = math.log(1.+9.*float(x)/float(nTicks),nBase)
            yVal = padY + h*(1. - yVal)
            self.IDstore.append(self.c.create_line(padX-3,yVal,padX+3,yVal,fill="black",tags="vol_hist"))
            self.IDstore.append(self.c.create_text(padX-16,yVal,text="%.2f"%(float(x)/float(nTicks)),tags="vol_hist",font=("helvetica", 9, "normal")))
    
        # lines to connect the data points
        for x in range(self.nBins):
            origX1, origY1 = self.binToX(x),   math.log(1.+9.*float(self.hist[x]-m)/(M-m),10.)
            origX2, origY2 = self.binToX(x+1), math.log(1.+9.*float(self.hist[x+1]-m)/(M-m),10.)
            x1, x2 = padX + origX1 * w, padX + origX2 * w
            y1, y2 = padY+h*(1-origY1), padY+h*(1-origY2)
            self.IDstore.append(self.c.create_line((x1,y1,x2,y2), fill="#3333ff",tags="vol_hist"))

        return self.c

if __name__=="__main__":
    root = Tk()

    # generate some fake data
    data = [ gauss(10,100) for x in range(500) ]

    v = VolumeHist(data,nBins=16,root=root,W=600,H=300,padX=30,padY=30)

    cc = v.toCanvas()

    root.mainloop()


#TODO:
#
#[2:27:18 PM] David Giesen: Yeah, my only suggestion on that X axis would be to go down to 2 sig figs
# (and maybe red for negative and black for +'ve?) and have a tooltip that pops up with a more accurate
# value (assuming a more accurate value is even of use...)
