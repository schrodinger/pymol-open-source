import os


# X-windows/TclTk specific paste hack

def get_selection():
	ftmp = os.environ['HOME']+"/.pymol-tmp"
	cmd = os.environ['PYMOL_PATH']+"/modules/paste.com "+ ftmp
	result = os.system(cmd)
	f=open(ftmp)
	sele = f.readlines()
	f.close()
	return sele

print get_selection()

