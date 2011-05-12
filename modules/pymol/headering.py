import os
import re
import sys
import string
import struct
import pymol
import math
from xray import space_group_map


if sys.byteorder == "big":
    MACH_ARCH_CODE_INT = 1
    MACH_ARCH_CODE_FLT = 1
else:
    MACH_ARCH_CODE_INT = 4
    MACH_ARCH_CODE_FLT = 4

class baseHeader:
    """
    A simple base class used to parse biological data headers
    """
    def __init__(self,filename):
        self.filename  = filename
        self.cols      = []
        self.byteorder_int = None
        self.byteorder_flt = None
        self.wordsize  = None

    def parseFile(self):
        pass
    def getColumns(self):
        return self.cols
    def getColumnsOfType(self,targetType):
        pass
    def checkFile(self):
        return os.path.isfile(self.filename)

    def guessCols(self,mapType):
        fC = self.getColumnsOfType("F")
        pC = self.getColumnsOfType("P")

        looksLike, FCol, PCol = [None]*3

        for prg in ("refmac", "phenix", "phenix_no_fill", "buster"):
            #
            # default (refmac,phenix,phenix_no_fill,buster) names are in the format:
            #  "2FoFc-amplitude 2FoFc-phase FoFc-ampl FoFc-phase"
            #
            defaultNames = pymol.cmd.get_setting_text("default_%s_names" % (prg))

            if len(defaultNames):
                defaultNames = string.split(defaultNames)
                if len(defaultNames):
                    if mapType=="2FoFc":
                        defaultF, defaultP = defaultNames[0], defaultNames[1]
                    elif mapType=="FoFc":
                        defaultF, defaultP = defaultNames[2], defaultNames[3]
                else:
                    print "Error: Please provide the setting 'default_%s_names' a comma separated string" % (prg)
                    print "       with the values for 2FoFc and FoFc, amplitude and phase names, respectively."
                    return [None]*3
            else:
                print "Error: Please provide the setting 'default_%s_names' a comma separated string" % (prg)
                print "       with the values for 2FoFc and FoFc, amplitude and phase names, respectively."
                return [None]*3
                        

            # defaultF =~ "FWT", "2FOFCFWT", "2FOFCWT_no_fill", "2FOFCWT", "DELFWT", "FOFCWT", ..
            for curFCol in fC:
                if curFCol.endswith(defaultF):
                    # found F, now look for matching P
                    pfx = string.rstrip(curFCol,defaultF)
                    curPCol = pfx+defaultP

                    ## print "curFCol = %s" % curFCol
                    ## print "curPCol = %s" % curPCol
                    
                    if curPCol in pC:
                        # found perfectly matching columns
                        looksLike = prg
                        FCol = curFCol
                        PCol = curPCol
                        return (FCol, PCol, looksLike)
        return [None]*3
                            
        

class CNSHeader(baseHeader):
    """
    CNS format
    """
    def __init__(self,filename):
        ## print "CNSHeader cstr"
        baseHeader.__init__(self,filename)
        self.parseFile()
    
class CIFHeader(baseHeader):
    """
    mmCIF
    """
    def __init__(self,filename):
        baseHeader.__init__(self,filename)
        self.parseFile()

    def parseFile(self):
        if self.checkFile():
            try:
                inFile = open(self.filename,'rb')
                data = inFile.readlines()
                inFile.close()

                in_loop = False
                
                curLine = data.pop(0)
                while curLine!=None:
                    if in_loop:
                        if curLine.startswith("_refln."):
                            self.cols.append(string.split(curLine,".",1)[1].strip())
                        else:
                            in_loop=False
                    else:
                        if curLine.startswith("loop_"):
                            in_loop=True
                    if len(data):
                        curLine = data.pop(0)
                    else:
                        curLine = None


            except IOError, e:
                print "Error-CIFReader: Couldn't read '%s' for input." % (self.filename)

class MTZHeader(baseHeader):
    HEADER_KEYWORDS = {
	"VERS"    : "VERS",
	"TITLE"   : "TITLE",
	"NCOL"    : "NCOL",
	"CELL"    : "CELL",
	"SORT"    : "SORT",
	"SYMINF"  : "SYMINF",
	"SYMM"    : "SYMM",
	"RESO"    : "RESO",
	"VALM"    : "VALM",
        "NDIF"    : "NDIF",
	"COL"     : "COL",
	"PROJECT" : "PROJECT",
	"CRYSTAL" : "CRYSTAL",
	"DATASET" : "DATASET",
	"DCELL"   : "DCELL",
	"DWAVEL"  : "DWAVEL",
	"BATCH"   : "BATCH",
        }

    """
    MTZ/CCP4
    """
    def __init__(self,filename):
        baseHeader.__init__(self,filename)

        self.version   = None
        self.title     = None
        self.ncol      = None
        self.nrefl     = None
        self.nbatches  = None
        self.cell      = None
        self.sort      = None
        self.syminf    = None
        self.symm      = []
        self.reso_min  = None
        self.reso_max  = None
        self.valm      = None
        self.ndif      = None
        self.col       = None
        self.project   = None
        self.crystal   = None
        self.dataset   = None
        self.dcell     = None
        self.dwavel    = None
        self.batch     = None

        self.datasets  = {}

        self.parseFile()
        self.format_cols()

    def format_cols(self,colType=None):
        """
        updates self.cols to a list of cols of a
        given MTZ type
        """
        self.cols = []

        # UNIQUE:
        # crystal_name/dataset_name/col_label"
        s = "/"
        c = []
        for key in self.datasets:
            c=[]
            ## print "KEY = %s" % key
            # user wants all columns
            if colType==None:
                c = self.datasets[key]["cols"].keys()
            else:
                # user wants a specfic column
                for tmpCol in self.datasets[key]["cols"].keys():
                    if self.datasets[key]["cols"][tmpCol]["type"]==colType:
                        c.append(tmpCol)
            ## print "Columns of type %s are: " % colType
            ## print c
            ## print "Now formatting...."
            if "crystal" in self.datasets[key]:
                curCryst = [ self.datasets[key]["crystal"] ] * len(c)
            else:
                curCryst = "nameless_crystal" * len(c)
            ## print "CurCrystal is: %s." % curCryst
            datName  = [ self.datasets[key]["name"] ] * len(c)
            self.cols.extend(
                map(lambda x: "/".join(x), 
                    zip(curCryst,datName,c)))

    def getColumns(self):
        self.format_cols()
        return baseHeader.getColumns(self)

    def getColumnsOfType(self,colType):
        self.format_cols(colType)
        return baseHeader.getColumns(self)
        
    def get_byteorder(self):
        f = open(self.filename, 'rb')
        f.seek(8)
        d = struct.unpack("BBBB",f.read(4))
        f.close()
        return (d[1]>>4) & 0x0f, (d[0]>>4) & 0x0f

    def authorStamp(self,f):
        # seek the machine stamp
        f.seek(8,0)
        # read machine stamp
        (B0,B1,B2,B3) = struct.unpack("BBBB", f.read(4))

        # determines big or little endian:
        # double, float, int, char
        d = (B0 & 0xF0) >> 4
        f = (B0 & 0x0F)
        i = (B1 & 0xF0) >> 4
        c = (B1 & 0x0F)

        if d==1:
            # big endian
            self.byteorder_flt = ">"
        elif d==4:
            # little endian
            self.byteorder_flt = "<"

        if i==1:
            # big
            self.byteorder_int = ">"
            self.wordsize = 1
        elif i==4:
            self.byteorder_int = "<"
            self.wordsize = 4
        
    def parseFile(self):
        if self.checkFile():
            f = open(self.filename,'rb')

            # sets authoring machine's details
            self.authorStamp(f)

            # get file size
            f.seek(0,2)
            file_len = f.tell()
            # get header offset
            f.seek(4)
            if self.wordsize!=None:
                (header_start,) = struct.unpack(self.byteorder_int+"i", f.read(4))
            else:
                print "Warning: Byte order of file unknown.  Guessing header location."
                (header_start,) = struct.unpack("i", f.read(4))

            # bAdjust is the byte adjustment to compensate for
            # older machines with smaller atomic sizes
            host_word_size = len(struct.pack('i',0))
            author_word_size = self.wordsize
            if host_word_size<author_word_size:
                bAdjust = author_word_size / host_word_size
            elif host_word_size>author_word_size:
                bAdjust = author_word_size * host_word_size
            else:
                bAdjust = host_word_size
                
            header_start  = (header_start-1) * (bAdjust)

            if file_len<header_start:
                print "Error: File '%s' cannot be parsed because PyMOL cannot find the header.  If you think"
                print "       PyMOL should be able to read this, plese send the file and this mesage to "
                print "       help@schrodinger.com.  Thanks!"

            f.seek(header_start)

            curLine = struct.unpack("80s", f.read(80))[0]
            
            while not (curLine.startswith("END")):
                # yank field identifier
                (field, tokens) = string.split(curLine," ",1)
                
                H = MTZHeader.HEADER_KEYWORDS
                try:
                    if field.startswith(H["VERS"]):
                        self.version = tokens
                    elif field.startswith(H["TITLE"]):
                        self.title = tokens
                    elif field.startswith(H["NCOL"]):
                        (self.ncols, self.nrefl, self.nbatches) = string.split(tokens)
                    elif field.startswith(H["CELL"]):
                        tokens = string.split(tokens)
                        self.cell_dim    = tokens[:3]
                        self.cell_angles = tokens[3:]
                    elif field.startswith(H["SORT"]):
                        self.sort = string.split(tokens)
                    elif field.startswith(H["SYMINF"]):
                        tokens = string.split(tokens,maxsplit=4)
                        self.nsymmop  = tokens[0]
                        self.nprimop  = tokens[1]
                        self.lattice  = tokens[2]
                        self.groupnum = tokens[3]

                        # try official MTZ format
                        # first 10 chars are the space_group
                        # from 10 on (6-chars) are the pt group
                        self.space_group = string.strip(tokens[4][:10],"'")
                        self.pt_group = string.strip(string.strip(tokens[4][10:],"'"))

                        # valid format for MTZ space group?
                        if self.space_group not in space_group_map.values():
                            # invalid format; try to guess
                            # SPCGP => S P C GP; maybe in the keys?
                            if self.space_group in space_group_map:
                                self.space_group = space_group_map[self.space_group]

                        # new refmac; smarter parsing
                        if self.space_group not in space_group_map.values():
                            m = re.match("'(?P<spc_gp>.*)' +(?P<pt_gp>.*)", tokens[4])
                            self.space_group, self.pt_group = m.group('spc_gp'), m.group('pt_gp')
                            if self.space_group in space_group_map:
                                self.space_group = space_group_map[self.space_group]

                        # still?
                        try:
                            if self.space_group not in space_group_map.values():
                                # eek, even worse (Buster); try to guess
                                (self.space_group, self.pt_group) = string.split(tokens[4])
                                self.space_group = string.strip(self.space_group,"'")
                                self.pt_group = string.strip(self.pt_group,"'")

                                if self.space_group in space_group_map:
                                    self.space_group = space_group_map[self.space_group]
                        except Exception as e:
                            pass
                        
                    elif field.startswith(H["SYMM"]):
                        tokens = string.split(tokens)
                        tokens = map(lambda x: string.strip(string.strip(x, ",")),tokens)
                        self.symm.append(tokens)
                    elif field.startswith(H["RESO"]):
                        self.reso_min, self.reso_max = string.split(tokens)
                        self.reso_min = math.sqrt(1./float(self.reso_min))
                        self.reso_max = math.sqrt(1./float(self.reso_max))
                    elif field.startswith(H["VALM"]):
                        self.missing_flag = tokens
                    elif field.startswith(H["COL"]):
                        (lab,typ,m,M,i) = string.split(tokens)
                        if i not in self.datasets:
                            self.datasets[i] = {}
                            self.datasets[i]["cols"]  = {}
                            
                        self.datasets[i]["cols"][lab] = {}
                        self.datasets[i]["cols"][lab]["type"] = typ
                        self.datasets[i]["cols"][lab]["min"]  = m
                        self.datasets[i]["cols"][lab]["max"]  = M
                        self.datasets[i]["cols"][lab]["id"]   = i
                    elif field.startswith(H["NDIF"]):
                        self.ndif = int(tokens)
                    elif field.startswith(H["PROJECT"]):
                        (i,proj) = string.split(tokens,maxsplit=1)
                        if i not in self.datasets:
                            self.datasets[i] = {}
                        self.datasets[i]["project"] = string.strip(proj)
                    elif field.startswith(H["CRYSTAL"]):
                        # can have multiple crystals per dataset?
                        # if, so not supported (overwritten) here.
                        (i,cryst) = string.split(tokens,maxsplit=1) 
                        if i not in self.datasets:
                            self.datasets[i] = {}
                        self.datasets[i]["crystal"] = string.strip(cryst)
                    elif field.startswith(H["DATASET"]):
                        (i,d) = string.split(tokens,maxsplit=1)
                        if i not in self.datasets:
                            self.datasets[i] = {}
                        self.datasets[i]["name"] = string.strip(d)
                    elif field.startswith(H["DCELL"]):
                        (i,x,y,z,a,b,g) = string.split(tokens)
                        if i not in self.datasets:
                            self.datasets[i] = {}
                        self.datasets[i]['x'] = x
                        self.datasets[i]['y'] = y
                        self.datasets[i]['z'] = z
                        self.datasets[i]['alpha'] = a
                        self.datasets[i]['beta'] = b
                        self.datasets[i]['gamma'] = g
                    elif field.startswith(H["DWAVEL"]):
                        (i,wl) = string.split(tokens)
                        if i not in self.datasets:
                            self.datasets[i] = {}
                        self.datasets[i]["wavelength"] = wl
                    elif field.startswith(H["BATCH"]):
                        self.batch = tokens
                    else:
                        print "Error Parsing MTZ Header: bad column name: '%s'" % field
                except ValueError:
                    print "Error: Parsing MTZ Header poorly formatted MTZ file"
                    print "       bad field: '%s'" % field
                    pass

                curLine = struct.unpack("80s", f.read(80))[0]


if __name__=="__main__":

    c = CIFHeader("test.cif")
    print c.getColumns()

    m = MTZHeader("test.mtz")
    print m.getColumns()
