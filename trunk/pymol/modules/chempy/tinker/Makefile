all: parm_simple parm99_simple


parm_simple: make_simple.py
	python make_simple.py > parm_simple.dat

parm99_simple: parm_simple.dat
	python make99_simple.py 

