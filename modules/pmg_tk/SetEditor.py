#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2006 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-* Kenneth Lind
#Z* -------------------------------------------------------------------

# Master editor for all pymol settings
# (includes filter)

from Tkinter import *

# nItem defines the number of label/entry widgets displayed on screen
# for scrolling.  Change to smaller number if performance is poor.
nItem = 15

class SetEditor:

    def __init__(self, app):
        self.pymol = app.pymol
        self.cmd = app.pymol.cmd
        
        #====================
        # get a list of settings and create dictionary with setting:value items
        self.list = self.pymol.setting.get_name_list()
        self.list.sort()
        self.index = {}
        for i in self.list:
            self.index[i] = self.cmd.get(i)

        #====================
        # create main frames and scale widget to use as scrollbar
        top = Toplevel()
        top.title( "PyMOL Settings" )
        top.geometry( "+100+200" )
        l = Label(top, text="Edit item and hit <Enter> to set value", bd=1, relief=RAISED )
        l.pack( side=TOP, fill=X, expand=1 )

        f1 = Frame(top)
        f1.pack( side=TOP, fill=BOTH, expand=1 )
        frame = Frame(f1)
        frame.pack( side=LEFT, fill=Y, expand=1 )

        self.scale = Scale( f1, to=len(self.list)-nItem, showvalue=0, command=self.updateLabels )
        self.scale.pack( side=RIGHT, fill=Y )

        #====================
        # create nItem label and entry widgets to hold setting name and values
        # the text in these widgets will change as the scaler is moved, 
        # giving a scrollbar effect
        self.labels = []
        self.values = []
        for i in range(nItem):
            l = Label( frame, text=self.list[i], width=30, anchor=E )
            l.grid( row=i, column=0, sticky=E )
            self.labels.append(l)

            e = Entry( frame, width=30 )
            e.insert( END, self.index[ self.list[i] ] )
            e.grid( row=i, column=1, sticky=W+E )
            e.bind( "<Return>", lambda event, i=i, e=e, s=self: s.onSet( i, e ) )
            self.values.append(e)

        #====================
        # create filter frame
        f2 = Frame( top, bd=1, relief=SUNKEN )
        f2.pack( side=BOTTOM, fill=X, expand=1 )
        l = Label( f2, text="Filter:", width=5, anchor=E )
        l.grid ( row=0, column=0 )
        self.filter = Entry( f2, width=48 )
        self.filter.bind( "<KeyRelease>", self.onFilter )
        self.filter.grid (row=0, column=1, sticky=W+E )
        b = Button( f2, text="Filter", width=5, command=self.onFilter )
        b.grid (row=0, column=2 )
        b = Button( f2, text="Reset", width=5, command=self.onResetFilter )
        b.grid (row=0, column=3 )

        top.focus_set()

    def updateLabels(self, *event):
        """update the labels and values in the window whenever the
        scale widget is 'scrolled', or whenever the filter is applied
        """
        if len(self.list) < nItem:
            self.scale.set(0)

        pos = self.scale.get()

        for i in range(nItem):
            try:
                self.labels[i].config( text=self.list[i+pos] )
            except:
                self.labels[i].config( text="" )

            self.values[i].delete(0,END)
            try:
                self.values[i].insert( END, self.index[ self.list[i+pos] ] )
            except:
                pass

    def onSet(self, offset, entry):
        """set the pymol setting to value in the entry widget where
        the enter key was pressed
        """
        try:
            lab = self.list[ self.scale.get()+offset ]
        except:
            return      # if trying to edit a blank entry

        val = entry.get()
        origVal = self.cmd.get( lab )
        try:
            self.cmd.set( lab, val, quiet=0, log=1 )
            self.index[lab] = val
        except:
            entry.delete( 0,END )
            entry.insert( END, origVal )

    def onFilter(self, *event):
        """get list of pymol settings that match filter"""
        val = self.filter.get()
        if not val: self.onResetFilter() # WLD
        newL = []
        for l in self.pymol.setting.get_name_list():
            if l.find( val ) != -1:
                newL.append( l )

        self.list = newL
        self.scale.set(0)
        self.scale.config( to=len(self.list)-nItem )
        self.updateLabels()

    def onResetFilter(self):
        """reset to full list of pymol settings"""
        self.list = self.pymol.setting.get_name_list()
        self.list.sort()
        self.scale.set(0)
        self.scale.config( to=len(self.list)-nItem )
        self.updateLabels()

