
import Tkinter as Tk

class ModalWindow(Tk.Toplevel):
    # 
    # Makes a toplevel window that is model (won't let you access the parent window when this window is open)
    # self.pack_and_display() = packs any buttons in the window and displays the window
    #
    # parent = window to display on top of
    # buttons = True to include a default set of "Accept" and "Cancel" buttons in the window
    # user_accept = the function that the "Accept" button should call before destroying itself
    # user_cancel = the function that the "Cancel" button should call before destroying itself
    #
    def __init__(self, parent, title='modal window', buttons=False, user_accept=False, user_cancel=False, *posargs, **kwargs):
        Tk.Toplevel.__init__(self, *posargs, **kwargs)
        self.parent_window = parent
        self.buttons = buttons
        self.user_accept = user_accept
        self.user_cancel = user_cancel
        self.parent_window = parent
        self.title(title)
        self.wm_protocol('WM_DELETE_WINDOW', self.cancel)
        self.withdraw()

    def cancel(self):
        # The user should define a function and assign it to the modal window as follows:
        # bob = modal_window(parent, user_cancel=function)
        if self.user_cancel:
            self.user_cancel()
        self.destroy()

    def accept(self,event):
        # The user should define a function and assign it to the modal window as follows:
        # bob = modal_window(parent, user_accept=function)
        if self.user_accept:
            self.user_accept(event)
        self.grab_release()
        self.destroy()

    def display(self, grab=True):
        #print "ModalWindow::display"
        # Don't include the height and width in the geometry string or else the window doesn't automatically resize
        self.update_idletasks()
        self.deiconify()
        # Transient makes this window not show up in the taskbar, and iconify, etc. with parent
        self.transient(self.parent_window)
        # focus_set makes this window active, which Windows doesn't always do without it
        self.focus_set()

        # The grab_set/wait_window combo means that this window stays on top and the parent
        # window doesn't respond to any mouseclicks
        if grab:
            # We may want to set grab=False if this window has widgets that generate another toplevel window 
            # (such as a DropDownCombobox)
            self.grab_set()
            #self.parent_window.wait_visibility(self.parent_window)
        self.parent_window.wait_window(self)
    
    def pack(self,x,y):
        if self.buttons:
            frame = Tk.Frame(self)
            accept_button = Tk.Button(frame, text='Accept', command=self.accept)
            cancel_button = Tk.Button(frame, text='Cancel', command=self.cancel)
            frame.pack(fill=Tk.X)
            accept_button.pack(side=Tk.LEFT)
            cancel_button.pack(side=Tk.RIGHT)
            self.buttons = False
        #print "Pack=(%d,%d)" % (x,y)

    def pack_and_display(self,event=None,grab=True):
        #print "ModalWindow::pack_and_display"
        wgt = event.widget
        (x,y) = wgt.canvasx(event.x), wgt.canvasy(event.y)
        self.pack(x,y)
        self.display(grab=grab)

def get_val():
    myvar.set(entry.get())


if __name__=="__main__":
    root = Tk.Tk()
    myvar = Tk.StringVar()
    myvar.set('Hello')
    label = Tk.Label(textvariable=myvar)
    mywin = ModalWindow(root, buttons=True, user_accept=get_val)
    button = Tk.Button(root, text='Show window', command=mywin.pack_and_display)
    entry = Tk.Entry(mywin)
    entry.pack()
    label.pack()
    button.pack()
    root.mainloop()

