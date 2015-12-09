from __future__ import print_function

try:
    import cStringIO
except ImportError:
    import io as cStringIO

import sys, os
import time
import re
from os.path import dirname
from subprocess import Popen, PIPE
from distutils import dir_util

def create_all(generated_dir, pymoldir="."):
    '''
    Generate various stuff
    '''
    create_shadertext(
            os.path.join(pymoldir, "data", "shaders"),
            generated_dir,
            "shadertext.txt",
            os.path.join(generated_dir, "ShaderText.h"),
            os.path.join(generated_dir, "ShaderText.c"))
    create_buildinfo(generated_dir, pymoldir)

class openw(object):
    """
    File-like object for writing files. File is actually only
    written if the content changed.
    """
    def __init__(self, filename):
        if os.path.exists(filename):
            self.out = cStringIO.StringIO()
            self.filename = filename
        else:
            dir_util.mkpath(os.path.dirname(filename))
            self.out = open(filename, "w")
            self.filename = None
    def close(self):
        if self.out.closed:
            return
        if self.filename:
            oldcontents = open(self.filename).read()
            newcontents = self.out.getvalue()
            if oldcontents != newcontents:
                self.out = open(self.filename, "w")
                self.out.write(newcontents)
        self.out.close()
    def __getattr__(self, name):
        return getattr(self.out, name)
    def __enter__(self):
        return self
    def __exit__(self, *a, **k):
        self.close()
    def __del__(self):
        self.close()

def create_shadertext(shaderdir, shaderdir2, inputfile, outputheader, outputfile):

    outputheader = openw(outputheader)
    outputfile = openw(outputfile)

    with open(os.path.join(shaderdir, inputfile)) as f:
        for l in f:
            lspl = l.split()
            if len(lspl)==0:
                continue
            if lspl[0] == "read":
                if len(lspl)!=3:
                    outputfile.write("/* WARNING: read doesn't have variable name and file name argument lspl=%s */\n" % lspl)
                else:
                    varname = lspl[1]
                    filename = lspl[2]
                    outputheader.write("extern const char* %s;\n" % varname)
                    outputfile.write("const char* %s =\n" % varname)
                    sd = shaderdir
                    if not os.path.exists(os.path.join(shaderdir, filename)):
                        sd = shaderdir2
                    with open(os.path.join(sd, filename)) as f2:
                        for l2 in f2:
                            st = l2.strip("\n")
                            if len(st)>0:
                                #if st[0] != '#':
                                outputfile.write("\"%s\\n\"\n" % st.replace('"', r'\"'))
                        outputfile.write(";\n") # end of variable definition
            else:
                outputheader.write("%s\n" % l.strip())
                outputfile.write("%s\n" % l.strip())

    outputheader.close()
    outputfile.close()

def create_buildinfo(outputdir, pymoldir='.'):

    try:
        sha = Popen(['git', 'rev-parse', 'HEAD'], cwd=pymoldir,
                stdout=PIPE).stdout.read().strip().decode()
    except OSError:
        sha = ''

    rev = 0
    try:
        for line in Popen(['svn', 'info'], cwd=pymoldir, stdout=PIPE).stdout:
            if line.startswith(b'Last Changed Rev'):
                rev = int(line.split()[3])
    except OSError:
        pass

    with openw(os.path.join(outputdir, 'PyMOLBuildInfo.h')) as out:
        print('''
#define _PyMOL_BUILD_DATE %d
#define _PYMOL_BUILD_GIT_SHA "%s"
#define _PyMOL_BUILD_SVN_REV %d
        ''' % (time.time(), sha, rev), file=out)

if __name__ == "__main__":
    create_shadertext(*sys.argv[1:6])
    create_buildinfo(dirname(sys.argv[4]), dirname(dirname(sys.argv[1])))
