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

# pm.py 
# Python interface module for PyMol
#
# **This is the only module which should be/need be imported by 
# **PyMol Programs

import re
import pmx
import string
import traceback

def export_dots(a,b):
	return pmx.export_dots(a,int(b))

def turn(a,b):
	pmx.turn(a,b)

def render():
	pmx.render()

def ray():
	pmx.render()

def zoom(a):
	pmx.zoom(a)

def frame(a):
	pmx.frame(int(a))

def move(a,b):
	pmx.move(a,b)

def clip(a,b):
	pmx.clip(a,b)

def origin(a):
	pmx.origin(a)

def refresh():
	pmx.refresh()

def set(a,b):
	pmx.set(a,b)

def turn(a,b):
	pmx.turn(a,b)

def reset():
	pmx.reset(0)
	
def delete(a):
	pmx.delete(a)


def quit():
	pmx.quit()

def png(a):
	pmx.png(a)

def mclear():
	pmx.mclear()

def mstop():
	pmx.mplay(0)

def mplay():
	pmx.mplay(1)

def mray():
	pmx.mplay(2)

def viewport(a,b):
	pmx.viewport(int(a),int(b))

def mdo(a,b):
	pmx.mdo(int(a)-1,b)

def run(*args):
# dummy. This functionality now provided directly by the parser.
	pass
			
def load(*args):
	ftype = 0
	if re.search("\.pdb$",args[0]):
		ftype = 0
	elif re.search("\.mol$",args[0]):
		ftype = 1
	if len(args)==1:
		oname = re.sub("[^/]*\/","",args[0])
		oname = re.sub("\.pdb|\.mol","",oname)
		pmx.load(oname,args[0],-1,ftype)
	elif len(args)==2:
		oname = string.strip(args[1])
		pmx.load(oname,args[0],-1,ftype)
	elif len(args)==3:
		oname = string.strip(args[1])
		pmx.load(oname,args[0],int(args[2])-1,ftype)
	else:
		print "argument error."

def read_mol(*args):
	ftype = 3
	if len(args)==2:
		oname = string.strip(args[1])
		pmx.load(oname,args[0],-1,ftype)
	elif len(args)==3:
		oname = string.strip(args[1])
		pmx.load(oname,args[0],int(args[2])-1,ftype)
	else:
		print "argument error."

def select(*args):
	if len(args)==1:
		sel_cnt = pmx.get("sel_counter") + 1.0
		pmx.set("sel_counter","%1.0f" % sel_cnt)
		sel_name = "sel%02.0f" % sel_cnt
		sel = args[0]
	else:
		sel_name = args[0]
		sel = args[1]
	pmx.select(sel_name,sel)

def color(*args):
	if len(args)==2:
		pmx.color(args[0],args[1],0)
	else:
		pmx.color(args[0],args[1],1)	

def mpng(*args):
	if len(args)==1:
		pmx.mpng(args[0],1)	
	elif args[1]=='purge':
		pmx.mpng(args[0],0);
	else:
		pmx.mpng(args[0],1);

def show(*args):
	if len(args)==2:
		if repres.has_key(args[0]):		
			repn = repres[args[0]];
			pmx.showhide(args[1],repn,1);
	elif args[0]=='all':
			pmx.showhide("!",0,1);

def hide(*args):
	if len(args)==2:
		if repres.has_key(args[0]):		
			repn = repres[args[0]];
			pmx.showhide(args[1],repn,0);
	elif args[0]=='all':
			pmx.showhide("!",0,0);

def mmatrix(a):
	if a=="clear":
		pmx.mmatrix(0)
	elif a=="store":
		pmx.mmatrix(1)
	elif a=="recall":
		pmx.mmatrix(2)

def enable(nam):
	pmx.onoff(nam,1);

def disable(nam):
	pmx.onoff(nam,0);

def mset(seq):
	output=[]
	input = string.split(seq," ")
	last = 0
	for x in input:
		if x[0]>"9" or x[0]<"0":
			if x[0]=="x":
				cnt = int(x[1:])-1
				while cnt>0:
					output.append(str(last))
					cnt=cnt-1
			elif x[0]=="-":
				dir=1
				cnt=last
				last = int(x[1:])-1
				if last<cnt:
					dir=-1
				while cnt!=last:
					cnt=cnt+dir
					output.append(str(cnt))
		else:
			val = int(x) - 1
			output.append(str(val))
			last=val
	pmx.mset(string.join(output," "))
	
keyword = { 
	'load'        : [load         , 1 , 3 , ',' , 0 ],
	'refresh'     : [refresh      , 0 , 0 , ',' , 0 ],
	'render'      : [render       , 0 , 0 , ',' , 0 ],
	'ray'         : [render       , 0 , 0 , ',' , 0 ],
	'select'      : [select       , 1 , 2 , '=' , 0 ],
	'set'         : [set          , 2 , 2 , '=' , 0 ],
	'frame'       : [frame        , 1 , 1 , ',' , 0 ],
	'turn'        : [turn         , 2 , 2 , ',' , 0 ],
	'move'        : [move         , 2 , 2 , ',' , 0 ],
	'clip'        : [clip         , 2 , 2 , ',' , 0 ],
	'show'        : [show         , 1 , 2 , ',' , 0 ],
	'hide'        : [hide         , 1 , 2 , ',' , 0 ],
	'disable'     : [disable      , 1 , 1 , ',' , 0 ],
	'enable'      : [enable       , 1 , 1 , ',' , 0 ],
	'delete'      : [delete       , 1 , 1 , ',' , 0 ],
	'zoom'        : [zoom         , 1 , 1 , ',' , 0 ],
	'origin'      : [origin       , 1 , 1 , ',' , 0 ],
	'color'       : [color        , 2 , 3 , ',' , 0 ],
	'quit'        : [quit         , 0 , 0 , ',' , 0 ],
	'mset'        : [mset         , 1 , 1 , ',' , 0 ],
	'png'         : [png          , 1 , 1 , ',' , 0 ],
	'mdo'         : [mdo          , 2 , 2 , ':' , 1 ],
	'mpng'        : [mpng         , 1 , 2 , ',' , 0 ],
	'mplay'       : [mplay        , 0 , 0 , ',' , 0 ],
	'mray'        : [mray         , 0 , 0 , ',' , 0 ],
	'mstop'       : [mstop        , 0 , 0 , ',' , 0 ],
	'mclear'      : [mclear       , 0 , 0 , ',' , 0 ],
	'mmatrix'     : [mmatrix      , 1 , 1 , ',' , 0 ],
	'viewport'    : [viewport     , 2 , 2 , ',' , 0 ],
	'reset'       : [reset        , 0 , 0 , ',' , 0 ],
	'export_dots' : [export_dots  , 2 , 2 , ',' , 0 ],
	'run'         : [run          , 1 , 2 , ',' , 2 ]
	}

repres = {
	'lines'       : 0,
	'sticks'      : 1,
	'dots'        : 2,
	'mesh'        : 3,
	'sphere'      : 4,
	'ribbon'      : 5
}



