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
#-* NOTE: Based on code by John E. Grayson which was in turn 
#-* based on code written by Doug Hellmann. 
#Z* -------------------------------------------------------------------

from Tkinter import *
from tkFileDialog import *

from AbstractApp import AbstractApp
from Setting import Setting

import Pmw
import sys, string
import pm
import re
import thread
import threading

class PMGApp(AbstractApp):

	appversion	  = '0.2'
	appname		 = 'PyMOL Molecular Graphics System'
	copyright	   = 'Copyright (C) 1998-2000 by Warren DeLano of\nDeLano Scientific. All rights reserved.'
	contactweb	  = 'http://www.pymol.org'
	contactemail	= 'warren@delanoscientific.com'

	def buttonAdd(self,frame,text,cmd):
		newBtn=self.createcomponent('button', (), None,
			Button,frame,text=text,highlightthickness=0,
			command=cmd,padx=0,pady=0)
		newBtn.pack(side=LEFT,fill=BOTH,expand=YES)
		
	def createButtons(self):
		row2 = self.createcomponent('row2', (), None,
			Frame,self.get_commandFrame(),bd=0)
		row2.pack(side=TOP,fill=BOTH,expand=YES)
		btn_reset = self.buttonAdd(row2,'Reset',pm.reset)
		btn_rtrace = self.buttonAdd(row2,'Ray Trace',lambda : threading.Thread(target=pm.ray,args=()).start())
		btn_reset = self.buttonAdd(row2,'Rock',pm.rock)

#		row3 = self.createcomponent('row3', (), None,
#			Frame,self.get_commandFrame(),bd=0)
#		row3.pack(side=TOP,fill=BOTH,expand=YES)
#		btn_reset = self.buttonAdd(row3,'Barf',self.quit)
#		btn_reset = self.buttonAdd(row3,'Now',self.quit)
#		btn_reset = self.buttonAdd(row3,'Eat',self.quit)
#		btn_reset = self.buttonAdd(row3,'Shrimp',self.quit)

		row1 = self.createcomponent('row1', (), None,
			Frame,self.get_commandFrame(),bd=0)
		row1.pack(side=TOP,fill=BOTH,expand=YES)
		btn_rewind = self.buttonAdd(row1,'|<',pm.rewind)
		btn_back = self.buttonAdd(row1,'<',pm.backward)
		btn_play = self.buttonAdd(row1,'Play',pm.mplay)
		btn_stop = self.buttonAdd(row1,'Stop',pm.mstop)
		btn_forward = self.buttonAdd(row1,'>',pm.forward)
		btn_last = self.buttonAdd(row1,'>|',pm.ending)
		btn_ccache = self.buttonAdd(row1,'Clear',pm.mclear)


	def createMain(self):
		self.command = StringVar()
		self.entry = self.createcomponent('entry', (), None,
									Entry,
									(self.get_dataArea(),),
									justify=LEFT,
									width=60,
									textvariable=self.command)
		self.entry.pack(side=BOTTOM,expand=NO,fill=X)
		self.entry.bind('<Return>',lambda event,w=self.command:
			(pm.do(w.get()),pm.dirty(),w.set('')))
		self.bind(self.entry, 'Command line entry field')

		self.output = self.createcomponent('output', (), None,
									Pmw.ScrolledText,
									(self.get_dataArea(),))
		self.output.after(50,self.update_feedback)
		self.output.pack(side=BOTTOM,expand=NO,fill=X)
		self.bind(self.entry, 'Output Area')

	def update_feedback(self):
		self.lineCount = self.lineCount + 1
		if self.lineCount > 10000:
			self.output.delete('0.0','%i.%i' % (self.lineCount-5000,0))
			self.lineCount=5000
		for a in pm.get_feedback():
			self.output.insert(END,"\n")
			self.output.insert(END,a)
			self.output.see(END)
			self.lineCount = self.lineCount + 1
		self.output.after(50,self.update_feedback)

	def createInterface(self):
			AbstractApp.createInterface(self)
			self.createButtons()
			self.createMain()
			self.lineCount = 0
			
	def quit_app(self):
		pm.quit()

	def file_open(self):
		ofile = askopenfilename(filetypes=[("PDB File","*.pdb"),("MOL File","*.mol")])
		if len(ofile):
			pm.load(ofile)

	def file_run(self):
		ofile = askopenfilename(filetypes=[("PyMOL Script","*.pml"),("Python Program","*.py")])
		if re.search("\.py$",ofile):
			pm.do("run "+ofile);		
		else:
			pm.do("@"+ofile);

	def file_save(self):
		sfile = asksaveasfilename(filetypes=[("PNG File","*.png")])
		if len(sfile):
			pm.png(sfile)

	def file_savemovie(self):
		sfile = asksaveasfilename(filetypes=[("Numbered PNG Files","*.png")])
		if len(sfile):
			pm.mpng(sfile)

	def createMenuBar(self):
		self.menuBar.addmenuitem('Help', 'command',
							 'Get information on application', 
							 label='About...', command=self.showAbout)

		self.toggleBalloonVar = IntVar()
		self.toggleBalloonVar.set(1)
		self.setting = Setting()
               
		self.menuBar.addmenuitem('Help', 'checkbutton',
							 	'Toggle balloon help',
							 	label='Balloon help',
								variable = self.toggleBalloonVar,
								command=self.toggleBalloon)

		self.menuBar.addmenuitem('File', 'command', 'Open structure file.',
								label='Open...',
								command=self.file_open)

#		self.menuBar.addmenuitem('File', 'command', 'Open sequential files.',
#								label='Open Sequence...',
#								command=self.file_open)

		self.menuBar.addmenuitem('File', 'separator', '')

		self.menuBar.addmenuitem('File', 'command', 'Run program or script.',
								label='Run...',
								command=self.file_run)

		self.menuBar.addmenuitem('File', 'separator', '')

		self.menuBar.addmenuitem('File', 'command', 'Save current image.',
								label='Save Image...',
								command=self.file_save)

		self.menuBar.addmenuitem('File', 'command', 'Save all frames.',
								label='Save Movie...',
								command=self.file_savemovie)

		self.menuBar.addmenuitem('File', 'separator', '')

		self.menuBar.addmenuitem('File', 'command', 'Quit PyMOL',
								label='Quit',
								command=pm.quit)


		self.menuBar.addmenu('Movie', 'Movie Control')

		self.menuBar.addmenuitem('Movie', 'checkbutton',
							 	'Photorealistic images.',
							 	label='Ray Trace Frames',
								variable = self.setting.ray_trace_frames,
								command = lambda s=self: s.setting.update('ray_trace_frames'))

		self.menuBar.addmenuitem('Movie', 'checkbutton',
							 	'Save images in memory.',
							 	label='Cache Frames',
								variable = self.setting.cache_frames,
								command = lambda s=self: s.setting.update('cache_frames'))

		self.menuBar.addmenu('Display', 'Display Control')
      
		self.menuBar.addmenuitem('Display', 'command', 'Stereo On',
                               label='Stereo On',
                               command = lambda: pm.stereo("on"))

		self.menuBar.addmenuitem('Display', 'command', 'Stereo Off',
                               label='Stereo Off',
                               command = lambda: pm.stereo("off"))

		self.menuBar.addmenu('Options', 'Configuration Control')

		self.menuBar.addmenuitem('Options', 'checkbutton',
							 	'Superimpose all molecular states.',
							 	label='Show All States',
								variable = self.setting.all_states,
								command = lambda s=self: s.setting.update('all_states'))

		self.menuBar.addmenuitem('Options', 'checkbutton',
							 	'Disable perspective.',
							 	label='Orthoscopic',
								variable = self.setting.ortho,
								command = lambda s=self: s.setting.update('ortho'))

		self.menuBar.addmenuitem('Options', 'checkbutton',
							 	'Smooth raytracing.',
							 	label='Antialias',
								variable = self.setting.antialias,
								command = lambda s=self: s.setting.update('antialias'))




