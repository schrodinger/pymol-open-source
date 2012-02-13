import sys, os

write = sys.stdout.write

def create_shadertext(shaderdir, inputfile, outputheader, outputfile):


    with open(os.path.join(shaderdir, inputfile)) as f:
        while 1:
            l = f.readline()
            if not l:
                break
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
                        while 1:
                            l2 = f2.readline()
                            if not l2:
                                outputfile.write(";\n") # end of variable definition
                                break
                            st = l2.strip("\n")
                            if len(st)>0:
                                outputfile.write("\"%s\\n\"\n" % st)
            else:
                outputheader.write("%s\n" % l.strip())
                outputfile.write("%s\n" % l.strip())

if __name__ == "__main__":
 
    shaderdir = sys.argv[1]
    inputfile = sys.argv[2]
    outputheader = open(sys.argv[3], 'w')
    outputfile = open(sys.argv[4], 'w')

    create_shader_text(shaderdir, inputfile, outputheader, outputfile)
