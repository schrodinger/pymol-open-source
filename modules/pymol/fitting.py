#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC. 
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

if __name__=='pymol.fitting':
	
	import cmd
	from cmd import _cmd,lock,unlock
	import selector
	import os
	import pymol
	import string
	
	from cmd import _cmd,lock,unlock,Shortcut, \
		  DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error


	def cealign(target, mobile, target_state=1, mobile_state=1, verbose=0, _self=cmd):
		'''
DESCRIPTION

	"cealign" aligns two proteins using the CE algorithm.

USAGE 

	cealign target, mobile [, target_state [, mobile_state [, verbose ]]]
]
NOTES
        The original algorithm has two settings (D0, and D1) for cutoffs and one window size.  We use their default settings (D0=3.0 Angstroms; D1=4.0 Angstroms; Window size=8) which were empircally determined through tests for speed and accuracy.
	Reference: Shindyalov IN, Bourne PE (1998) Protein structure alignment by incremental combinatorial extension (CE) of the optimal path. Protein Engineering 11(9) 739-747.

	
EXAMPLES

	cealign protA////CA, protB////CA

	# fetch two proteins and align them
	fetch 1rlw 1rsy
	cealign 1rlw, 1rsy

SEE ALSO

	align, pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur, super
		'''
		#########################################################################
		# CE specific defines.	Don't change these unless you know
		# what you're doing.  See the documentation.
		#########################################################################
		# WINDOW SIZE
		# make sure you set this variable in cealign.py, as well!
		winSize = 8
		# FOR AVERAGING
		winSum = (winSize-1)*(winSize-2) / 2;
		# max gap size
		gapMax = 30

		quiet = not(int(verbose))

		from pymol import stored
		# make the lists for holding coordinates
		# partial lists
		stored.sel1, stored.sel2 = [], []
		# full lists
		stored.mol1, stored.mol2 = [], []

		# now put the coordinates into a list

		# -- REMOVE ALPHA CARBONS
		mobile = mobile + " and n. CA"
		target = target + " and n. CA"
		# -- REMOVE ALPHA CARBONS

		r = DEFAULT_ERROR
		
		# handle PyMOL's macro /// notation
		mobile = selector.process(mobile)
		target = selector.process(target)
				
		cmd.iterate_state(mobile_state, mobile, "stored.sel2.append([x,y,z])")
		cmd.iterate_state(target_state, target, "stored.sel1.append([x,y,z])")
		
		# full molecule
		mol1 = cmd.identify(target,1)[0][0]
		mol2 = cmd.identify(mobile,1)[0][0]
		
		# put all atoms from MOL1 & MOL2 into stored.mol1; gets the VLA of coords
		cmd.iterate_state(target_state, mol1, "stored.mol1.append([x,y,z])")
		cmd.iterate_state(mobile_state, mol2, "stored.mol2.append([x,y,z])")
		
		if ( len(stored.mol1) == 0 ):
			print "ERROR: Your first selection was empty."
			return
		if ( len(stored.mol2) == 0 ):
			print "ERROR: Your second selection was empty."
			return

		try:
			_self.lock(_self)

			# call the C function
			r = _cmd.cealign( _self._COb, stored.sel1, stored.sel2 )

			(aliLen, RMSD, rotMat) = r
			if not quiet:
				import pprint
				print "Aligned %d residues." % (aliLen)
				print "RMSD of Alignment: %f" % (RMSD)
				print "TTT Matrix:"
				pprint.pprint(rotMat)
			else:
				print "RMSD %f over %i residues" % (float(RMSD), int(aliLen))

			cmd.transform_selection( mol2, rotMat, homogenous=0 );
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return ( {"alignment_length": aliLen, "RMSD" : RMSD, "rotation_matrix" : rotMat } )


	def alignto(sel1,quiet=1):
		"""
DESCRIPTION

	NOTE: This feature is experimental and unsuspported.

	"alignto" aligns all other loaded objects to the given selected object
	using the CEalign algorithm.

USAGE

	alignto target [, quiet ]

EXAMPLE

        # fetch some calmodulins
	fetch 1cll 1sra 1ggz 1k95, async=0
	alignto 1cll

SEE ALSO

        align, super, cealign, fit, rms, rms_cur, intra_fit
		"""
		for x in cmd.get_names("objects"):
			if not quiet:
        	        	print "Aligning %s to %s" % (x, sel1)
			cealign( sel1, x )
				


	def super(mobile, target, cutoff=2.0, cycles=5,
			  gap=-1.5, extend=-0.7, max_gap=50, object=None,
			  matrix="BLOSUM62", mobile_state=0, target_state=0, 
			  quiet=1, max_skip=0, transform=1, reset=0,
			  seq=0.0, radius=12.0, scale=17.0, base=0.65,
			  coord=0.0, expect=6.0, window=3, ante=-1.0,
			  _self=cmd):
		
		'''
DESCRIPTION

	NOTE: This feature is experimental and unsupported.
	
	"super" performs a residue-based pairwise alignment followed by a
	structural superposition, and then carries out zero or more cycles
	of refinement in order to reject outliers.

USAGE 

	super mobile, target [, object=name ]

NOTES

	By adjusting various parameters, the nature of the initial
	alignment can be modified to include or exclude various factors
	including sequence similarity, main chain path, secondary &
	tertiary structure, and current coordinates.

EXAMPLE

	super protA////CA, protB////CA, object=supeAB

SEE ALSO

	align, pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur
	'''
		r = DEFAULT_ERROR
		mobile = selector.process(mobile)
		target = selector.process(target)
		if object==None: object=''
		matrix = str(matrix)
		if string.lower(matrix)=='none':
			matrix=''
		if len(matrix):
			mfile = cmd.exp_path("$PYMOL_DATA/pymol/matrices/"+matrix)
		else:
			mfile = ''		  
		# delete existing alignment object (if asked to reset it)
		try:
			_self.lock(_self)
			
			r = _cmd.align(_self._COb,mobile,"("+target+")",float(cutoff),
						   int(cycles),float(gap),float(extend),int(max_gap),
						   str(object),str(mfile),
						   int(mobile_state)-1,int(target_state)-1,
						   int(quiet),int(max_skip),int(transform),
						   int(reset),float(seq),
						   float(radius),float(scale),float(base),
						   float(coord),float(expect),int(window),
						   float(ante))
			
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def align(mobile, target, cutoff=2.0, cycles=5, gap=-10.0,
			  extend=-0.5, max_gap=50, object=None,
			  matrix="BLOSUM62", mobile_state=0, target_state=0,
			  quiet=1, max_skip=0, transform=1, reset=0, _self=cmd):
		
		'''
DESCRIPTION

	"align" performs a sequence alignment followed by a structural
	superposition, and then carries out zero or more cycles of
	refinement in order to reject outliers.

USAGE 

	align mobile, target [, object=name ]

ARGUMENTS

	mobile = string: atom selection for mobile atoms

	target = string: atom selection for target atoms

	object = string: name of alignment object to create
	
NOTES

	If object is specified, then align will create an object which
	indicates paired atoms and supports visualization of the alignment
	in the sequence viewer.

EXAMPLE

	align protA////CA, protB////CA, object=alnAB

SEE ALSO

	super, pair_fit, fit, rms, rms_cur, intra_rms, intra_rms_cur
		'''
		r = DEFAULT_ERROR
		mobile = selector.process(mobile)
		target = selector.process(target)
		matrix = str(matrix)
		if string.lower(matrix)=='none':
			matrix=''
		if len(matrix):
			mfile = cmd.exp_path("$PYMOL_DATA/pymol/matrices/"+matrix)
		else:
			mfile = ''
		if object==None: object=''
		# delete existing alignment object (if asked to reset it)
		try:
			_self.lock(_self)
			r = _cmd.align(_self._COb,mobile,"("+target+")",
						   float(cutoff),int(cycles),float(gap),
						   float(extend),int(max_gap),str(object),str(mfile),
						   int(mobile_state)-1,int(target_state)-1,
						   int(quiet),int(max_skip),int(transform),int(reset),
						   -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0)
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def intra_fit(selection, state=1, quiet=1, mix=0, _self=cmd):
		'''
DESCRIPTION

	"intra_fit" fits all states of an object to an atom selection
	in the specified state.	 It returns the rms values to python
	as an array.

USAGE 

	intra_fit selection [, state]

ARGUMENTS

	selection = string: atoms to fit

	state = integer: target state

PYMOL API

	cmd.intra_fit( string selection, int state )

EXAMPLES

	intra_fit ( name ca )

PYTHON EXAMPLE

	from pymol import cmd
	rms = cmd.intra_fit("(name ca)",1)

SEE ALSO

	fit, rms, rms_cur, intra_rms, intra_rms_cur, pair_fit
		'''
		# preprocess selection
		selection = selector.process(selection)
		#	
		r = DEFAULT_ERROR
		state = int(state)
		mix = int(mix)
		try:
			_self.lock(_self)
			r = _cmd.intrafit(_self._COb,"("+str(selection)+")",int(state)-1,2,int(quiet),int(mix))
		finally:
			_self.unlock(r,_self)
		if r<0.0:
			r = DEFAULT_ERROR
		elif not quiet:
			st = 1
			for a in r:
				if a>=0.0:
					if mix:
						print " cmd.intra_fit: %5.3f in state %d vs mixed target"%(a,st)
					else:
						print " cmd.intra_fit: %5.3f in state %d vs state %d"%(a,st,state)
				st = st + 1
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def intra_rms(selection, state=0, quiet=1, _self=cmd):
		'''
DESCRIPTION

	"intra_rms" calculates rms fit values for all states of an object
	over an atom selection relative to the indicated state.
	Coordinates are left unchanged.	 The rms values are returned as a
	python array.

EXAMPLE

	from pymol import cmd
	rms = cmd.intra_rms("(name ca)",1)
	print rms

PYMOL API

	cmd.intra_rms(string selection, int state)

SEE ALSO

	fit, rms, rms_cur, intra_fit, intra_rms_cur, pair_fit
		'''
		# preprocess selection
		selection = selector.process(selection)
		#	
		r = DEFAULT_ERROR
		state = int(state)
		try:
			_self.lock(_self)
			r = _cmd.intrafit(_self._COb,"("+str(selection)+")",int(state)-1,1,int(quiet),int(0))
		finally:
			_self.unlock(r,_self)
		if r<0.0:
			r = DEFAULT_ERROR
		elif not quiet:
			st = 1
			for a in r:
				if a>=0.0:
					print " cmd.intra_rms: %5.3f in state %d vs state %d"%(a,st,state)
				st = st + 1
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def intra_rms_cur(selection, state=0, quiet=1, _self=cmd):
		'''
DESCRIPTION

	"intra_rms_cur" calculates rms values for all states of an object
	over an atom selection relative to the indicated state without
	performing any fitting.	 The rms values are returned
	as a python array.

PYMOL API

	cmd.intra_rms_cur( string selection, int state)

PYTHON EXAMPLE

	from pymol import cmd
	rms = cmd.intra_rms_cur("(name ca)",1)

SEE ALSO

	fit, rms, rms_cur, intra_fit, intra_rms, pair_fit
		'''
		# preprocess selection
		selection = selector.process(selection)
		#	
		r = DEFAULT_ERROR
		state = int(state)
		try:
			_self.lock(_self)
			r = _cmd.intrafit(_self._COb,"("+str(selection)+")",int(state)-1,0,int(quiet),int(0))
		finally:
			_self.unlock(r,_self)
		if r<0.0:
			r = DEFAULT_ERROR
		elif not quiet:
			st = 1
			for a in r:
				if a>=0.0:
					print " cmd.intra_rms_cur: %5.3f in state %d vs state %d"%(a,st,state)
				st = st + 1
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def fit(mobile, target, mobile_state=0, target_state=0,
			  quiet=1, matchmaker=0, cutoff=2.0, cycles=0, object=None, _self=cmd):
		'''
DESCRIPTION

	"fit" superimposes the model in the first selection on to the model
	in the second selection.
	
USAGE

	fit mobile, target

EXAMPLES

	fit protA, protB

NOTES

	Only matching atoms in both selections will be used for the fit.
	
	Since atoms are matched based on all of their identifiers
	(including segment and chain identifiers), this command is only
	helpful when comparing very similar structures.

SEE ALSO

	align, super, pair_fit, rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
		'''
		r = DEFAULT_ERROR	   
		a=str(mobile)
		b=str(target)
		# preprocess selections
		a = selector.process(a)
		b = selector.process(b)
		#
		if object==None: object=''
		if int(matchmaker)==0:
			sele1 = "((%s) in (%s))" % (str(a),str(b))
			sele2 = "((%s) in (%s))" % (str(b),str(a))
		else:
			sele1 = str(a)
			sele2 = str(b)
		try:
			_self.lock(_self)
			r = _cmd.fit(_self._COb,sele1,sele2,2,
							 int(mobile_state)-1,int(target_state)-1,
							 int(quiet),int(matchmaker),float(cutoff),
							 int(cycles),str(object))
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def rms(mobile, target, mobile_state=0, target_state=0, quiet=1,
			  matchmaker=0, cutoff=2.0, cycles=0, object=None, _self=cmd):
		'''
DESCRIPTION

	"rms" computes a RMS fit between two atom selections, but does not
	tranform the models after performing the fit.

USAGE

	rms (selection), (target-selection)

EXAMPLES

	fit ( mutant and name ca ), ( wildtype and name ca )

SEE ALSO

	fit, rms_cur, intra_fit, intra_rms, intra_rms_cur, pair_fit	  
		'''
		r = DEFAULT_ERROR	   
		a=str(mobile)
		b=str(target)
		# preprocess selections
		a = selector.process(a)
		b = selector.process(b)
		#
		if object==None: object=''		
		if int(matchmaker)==0:
			sele1 = "((%s) in (%s))" % (str(a),str(b))
			sele2 = "((%s) in (%s))" % (str(b),str(a))
		else:
			sele1 = str(a)
			sele2 = str(b)
		try:
			_self.lock(_self)	
			r = _cmd.fit(_self._COb,sele1,sele2,1,
							 int(mobile_state)-1,int(target_state)-1,
							 int(quiet),int(matchmaker),float(cutoff),
							 int(cycles),str(object))
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def rms_cur(mobile, target, mobile_state=0, target_state=0,
				quiet=1, matchmaker=0, cutoff=2.0, cycles=0,
				object=None, _self=cmd):
		
		'''
DESCRIPTION

	"rms_cur" computes the RMS difference between two atom
	selections without performing any fitting.

USAGE

	rms_cur (selection), (selection)

SEE ALSO

	fit, rms, intra_fit, intra_rms, intra_rms_cur, pair_fit	  
		'''
		r = DEFAULT_ERROR	   
		a=str(mobile)
		b=str(target)
		# preprocess selections
		a = selector.process(a)
		b = selector.process(b)
		#
		if object==None: object=''			  
		if int(matchmaker)==0:
			sele1 = "((%s) in (%s))" % (str(a),str(b))
			sele2 = "((%s) in (%s))" % (str(b),str(a))
		else:
			sele1 = str(a)
			sele2 = str(b)
		try:
			_self.lock(_self)
			r = _cmd.fit(_self._COb,sele1,sele2,0,
							 int(mobile_state)-1,int(target_state)-1,
							 int(quiet),int(matchmaker),float(cutoff),
							 int(cycles),str(object))
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r

	def pair_fit(*arg, **kw):
		'''
DESCRIPTION

	"pair_fit" fits matched sets of atom pairs between two objects.

USAGE

	pair_fit selection, selection, [ selection, selection [ ... ]]

EXAMPLES

	# superimpose protA residues 10-25 and 33-46 to protB residues 22-37 and 41-54:
	
	pair_fit protA/10-25+33-46/CA, protB/22-37+41-54/CA

	# superimpose ligA atoms C1, C2, and C4 to ligB atoms C8, C4, and C10, respectively:
	
	pair_fit ligA////C1, ligB////C8, ligA////C2, ligB////C4, ligA////C3, ligB////C10
	
NOTES

	So long as the atoms are stored in PyMOL with the same order
	internally, you can provide just two selections.  Otherwise, you
	may need to specify each pair of atoms separately, two by two, as
	additional arguments to pair_fit.
	
	Script files are usually recommended when using this command.

SEE ALSO

	fit, rms, rms_cur, intra_fit, intra_rms, intra_rms_cur
		'''
		_self = kw.get('_self',cmd)
		r = DEFAULT_ERROR	   
		new_arg = []
		for a in arg:
			new_arg.append(selector.process(a))
		try:
			_self.lock(_self)	
			r = _cmd.fit_pairs(_self._COb,new_arg)
		finally:
			_self.unlock(r,_self)
		if _self._raising(r,_self): raise pymol.CmdException		 
		return r






