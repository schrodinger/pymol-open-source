import sys, os
from distutils import dir_util

def openw(filename):
    if not isinstance(filename, basestring):
        return filename
    dir_util.mkpath(os.path.dirname(filename))
    return open(filename, 'w')

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

if __name__ == "__main__":
    create_shader_text(*sys.argv[1:5])
