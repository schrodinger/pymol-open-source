import headering
import glob
import pprint
import sys
import os

if len(sys.argv)>1:
	p = sys.argv[1:]
else:
	print "Usage: pymol -cq ./headerTest.py -- MTZTestDirectory"
	sys.exit(1)


for curPath in p:
	print "Processing MTZ files in path: '%s'" % curPath

	l = glob.glob(curPath + os.sep + "*.mtz")

	d = {}

	for mtzFile in l:
	    print "Processing %s." % mtzFile
	    d[mtzFile] = headering.MTZHeader(mtzFile)

	for f in d:
	    print "Columns for file '%s':" % f
	    pprint.pprint(d[f].getColumns())
	    print ""

	for f in d:
	    print "Columns of type W ofr file '%s':" % f
	    pprint.pprint(d[f].getColumnsOfType("W"))
	    print ""
		
