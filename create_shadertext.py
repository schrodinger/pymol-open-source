import sys, os, cStringIO
from distutils import dir_util

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

def create_shadertext(shaderdir, inputfile, outputheader, outputfile):

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
                    with open(os.path.join(shaderdir, filename)) as f2:
                        for l2 in f2:
                            st = l2.strip("\n")
                            if len(st)>0:
                                outputfile.write("\"%s\\n\"\n" % st)
                        outputfile.write(";\n") # end of variable definition
            else:
                outputheader.write("%s\n" % l.strip())
                outputfile.write("%s\n" % l.strip())

    outputheader.close()
    outputfile.close()

if __name__ == "__main__":
    create_shadertext(*sys.argv[1:5])
