#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

# pmg.py 
# TKinter based gui for PyMol
# NOTE: must have treads support compiled into python to use this module
#
# **This is the only module which should be/need be imported by 
# **PyMol Programs

import pm
import thread
import sys

from Tkinter import *
import Pmw

def frame(root,side):
	w = Frame(root)
	w.pack(side=side, expand=YES,fill=BOTH)
	return w

def button(root, side, text, command = None):
	w = Button(root, text=text, command=command)
	w.pack(side=side, expand=YES, fill=BOTH)
	return w

class PMCalc(Frame):
	def __init__(self):
		root = Tk()
		Frame.__init__(self)
		self.pack(expand=YES, fill=BOTH)
		self.master.title('PyMOL GUI')
		self.master.iconname("pmgui")

		display = StringVar()
		entry = Entry(self,relief=SUNKEN,textvariable=display)
		entry.pack(side=TOP,expand=YES,fill=BOTH)
	
		for key in ("123","456","789","-0."):
			keyF=frame(self,TOP)
			for char in key:
				button(keyF,LEFT,char,lambda w=display,s=' %s '%char: w.set(w.get()+s))
	
		opsF = frame(self,TOP)
	
		for char in "+-*/=":
			if char == '=':
				btn = button(opsF,LEFT,char)
				btn.bind('<ButtonRelease-1>',
					lambda e, s=self, w=display: s.calc(w),'+')
			else:
				btn = button(opsF, LEFT, char,
				lambda w=display,c=char:w.set(w.get()+' '+c+' '))
		clearF = frame(self,BOTTOM)
		button(clearF,LEFT,'Clr', lambda w=display: w.set(''))

	def calc(self,display):
		try:
			display.set(`eval(display.get())`)
		except ValueError:
			display.set("ERROR")

class PMG:
	
	def __init__(self):
		self.root = Tk()
		pm.do("load ../test.mol")
		self.root.geometry('640x100+0+0')
		self.root.title("PyMOL Console")
		Pmw.initialise()		
		self.command = StringVar()
		self.entry = Entry(self.root,width=80,relief=SUNKEN,justify=LEFT,
				textvariable=self.command)
		self.entry.pack(side=BOTTOM,expand=NO,fill=BOTH)
		self.entry.bind('<Return>',lambda event,w=self.command:
			(pm.do(w.get()),pm.refresh(),w.set('')))
		self.root.mainloop()	

def main():
	sys.argv=["pymol"]
	PMG()
	
thread.start_new_thread(main,())




