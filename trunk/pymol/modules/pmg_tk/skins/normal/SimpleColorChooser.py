import os
import colorsys
from Tkinter import *
import ModalWindow

selectedColor = None

def setColor(event):
    global selectedColor
    w = event.widget
    X, Y = w.canvasx(event.x), w.canvasy(event.y)
    W,H = w.winfo_reqwidth(), w.winfo_reqheight()
    selectedColor = colorsys.hsv_to_rgb(X/W,Y/H,1.0)
    print "callback SimpleColorChooser::setColor-User selected a color, setting selectedColor to:\n\t", str(selectedColor)

class SimpleColorChooser(ModalWindow.ModalWindow):
    def pack(self,x,y):
        self.padX = self.padY = 20
        self.frame = Frame(self,width=360+self.padX,height=100+self.padY)
        self.canvas = Canvas(self.frame,width=360+self.padX,height=100+self.padY)
        self.canvas.pack()
        self.imgName = ""
        if "PYMOL_DATA" in os.environ:
            self.imgName = os.environ["PYMOL_DATA"] + "/pymol/hsv.ppm"
        else:
            self.imgName = "hsv.ppm"
        print "self.imgName = ", self.imgName
        self.img = PhotoImage(file=self.imgName)
        self.canvas.create_image((self.padX/2, self.padY/2),image=self.img,anchor=NW)
        self.canvas.bind("<Button-1>", self.accept)
        print "SimpleColorChooser::pack7"
        self.frame.pack()

            
if __name__=="__main__":

    root = Tk()
    a = SimpleColorChooser(root, title="Simple Color Chooser", user_accept=setColor)
    root.bind("<Button-1>", a.pack_and_display)
    root.wait_window(a)
    print "user selected", selectedColor
    root.mainloop()


#root = Tk()
#c=Canvas(root,width=360,height=100)
#c.pack()
#img = PhotoImage(file="hsv.ppm")
#c.create_image(0,0,image=img,anchor=NW)
