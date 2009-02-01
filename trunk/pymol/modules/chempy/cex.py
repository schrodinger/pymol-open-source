#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Scott Dixon, Metaphorics, LLC
#-* 
#-*
#Z* -------------------------------------------------------------------

"""
Author: Scott Dixon, Metaphorics, LLC
This source code is contributed to the public domain and may be freely
copied and distributed for research, profit, fun or any other reason,
with these restrictions: (1) unmodified or functionally equivalent code
derived from this code must contain this notice, (2) all derived code
must acknowledge the author and institution, and (3) no liability is
assumed by the author(s) for any use or misuse of this software.

CEX input routines.  Reads each CEX object into a test based tree.
Provides a CEX smiles interpreter class which can be specialized to create
appropriate molecule object """
import string

class CEXstream:
    """Input stream which read from file object"""
    (START, COMMENT, QUOTE, NOTQUOTE, GOTQUOTE, TAG, VALUE, END) = range(8)
    TAG_CHAR = string.letters + string.digits + "$_/"
    def __init__(self,file):
        self.file = file
        self.dt=None
        self.oldch = 0
        self.buff = ""
        self.p = 0
        self.len = 0
    def readEntry(self):
        """Read one tag<value> entry from stream"""
        # find nonblank character
        str = ""
        p = 0
        while 1:
            try:
                if self.buff[p] not in string.whitespace:
                    break
                p = p + 1
            except IndexError:
                self.buff = self.file.read(1000)
                p = 0
                if len(self.buff) == 0:
                    return (None, None)
        self.buff = self.buff[p:]
        if self.buff[0] == "|":
            self.buff = self.buff[1:]
            return ("|","")
        while 1:
            try:
                while 1:
                    p = string.index(self.buff,">") + 1
                    str = str + self.buff[:p]
                    self.buff = self.buff[p:]
                    if string.count(str,'"') %2 == 0:
                        break
            except (ValueError, IndexError):
                str = str + self.buff
                self.buff = self.file.read(1000)
                if len(self.buff)==0:
                    if string.find(str,"|") >= 0:
                        return ("|","")
                    else:
                        return (None, None)
            else: break    
        s = string.find(str,"<")
        if s < 0:
            return (None, None)
        else:
            return (str[:s],str[s+1:-1])
        
class CEXsmilesError(Exception):
    def __init__(self,smiles,p,msg):
        self.args="Smiles error: " + msg + "\n" + smiles + "\n" + p*" " + "^"
class CEXsmilesParser:
    """A simple CEX smiles parser adapted from Dave Weininger's C version in the
    CEX toolkit"""
    MX_NESTING=4096
    MX_RINGS=1000
    ptab = {"*":0, 
   "H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10,
   "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18, "K":19, "Ca":20,
   "Sc:":21, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30,
   "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39, "Zr":40,
   "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50,
   "Sb":51, "Te":52, "I":53, "Xe":54, "Cs":55, "Ba":56, "La":57, "Ce":58, "Pr":59, "Nd":60,
   "Pm":61, "Sm":62, "Eu":63, "Gd":64, "Tb":65, "Dy":66, "Ho":67, "Er":68, "Tm":69, "Yb":70,
   "Lu":71, "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79, "Hg":80,
   "Tl":81, "Pb":82, "Bi":83, "Po":84, "At":85, "Rn":86, "Fr":87, "Ra":88, "Ac":89, "Th":90,
   "Pa":91, "U":92, "Np":93, "Pu":94, "Am":95, "Cm":96, "Bk":97, "Cf":98, "Es":99, "Fm":100,
   "Md":101, "No":102, "Lr":103, "Rf":104, "Ha":105}
    stab = {0:"*", 
   1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 10:"Ne",
   11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 18:"Ar", 19:"K", 20:"Ca",
   21:"Sc:", 22:"Ti", 23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 30:"Zn",
   31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 36:"Kr", 37:"Rb", 38:"Sr", 39:"Y", 40:"Zr",
   41:"Nb", 42:"Mo", 43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd", 49:"In", 50:"Sn",
   51:"Sb", 52:"Te", 53:"I", 54:"Xe", 55:"Cs", 56:"Ba", 57:"La", 58:"Ce", 59:"Pr", 60:"Nd",
   61:"Pm", 62:"Sm", 63:"Eu", 64:"Gd", 65:"Tb", 66:"Dy", 67:"Ho", 68:"Er", 69:"Tm", 70:"Yb",
   71:"Lu", 72:"Hf", 73:"Ta", 74:"W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au", 80:"Hg",
   81:"Tl", 82:"Pb", 83:"Bi", 84:"Po", 85:"At", 86:"Rn", 87:"Fr", 88:"Ra", 89:"Ac", 90:"Th",
   91:"Pa", 92:"U", 93:"Np", 94:"Pu", 95:"Am", 96:"Cm", 97:"Bk", 98:"Cf", 99:"Es", 100:"Fm",
   101:"Md", 102:"No", 103:"Lr", 104:"Rf", 105:"Ha"}
    def sym2num(self,sym):
        try:
            return CEXsmilesParser.ptab[sym]
        except KeyError:
            return -1
    def num2sym(self,num):
        try:
            return CEXsmilesParser.stab[num]
        except KeyError:
            return ""
    def needquote(self,atnum):
        if atnum in (0,5,6,7,8,9,15,16,17,35,53): return 0
        else: return 1
    def __init__(self):
        self.atomN = 0
    def MakeAtom(self, atnum):
        print "Atom %d, atomic number %d" % (self.atomN, atnum)
        self.atomN = self.atomN + 1
        return self.atomN-1
    def MakeBond(self, at1, at2, bo):
        print "Bond between %d and %d, order %d" % (at1, at2,bo)
    def SetHcount(self, atom, count):
        print "Explicit H count %d for atom %d" % (count, atom)
    def SetFormalCharge(self, atom, charge):
        print "Charge for atom %d is %d" % (atom, charge)
    def SetAtomicMass(self, atom, mass):
        print "Mass from atom %d is %d" % (atom, mass)
    def parse(self,smiles):
        self.smiles=smiles + 3*"\0"  # guard zone for illegal smiles
        self.__init__()
        self.ringat = [None]*CEXsmilesParser.MX_RINGS
        self.fromat = [None]*CEXsmilesParser.MX_RINGS
        self.ringbo = [0]*CEXsmilesParser.MX_NESTING
        self.molname = ""
        lev = 0
        atnum = -1
        imph = -1
        bo = 0
        charge = 0
        quoted = 0
        mass = 0
        # adapted from Dave Wieninger's code in the CEX toolkits
        p = 0
        while p < len(self.smiles):
            pp = p + 1
            ch = self.smiles[p]
            if ch == "(":
                self.fromat[lev + 1] = self.fromat[lev]
                lev = lev + 1
            elif ch == ")": lev = lev - 1
            elif ch == "[":
                if quoted:
                    # error, no closing ]
                    raise CEXsmilesError(smiles,p,"No closing ]")
                else:
                    quoted = 1
                    if self.smiles[pp] in string.digits:
                        p = pp
                        while self.smiles[p+1] in string.digits:
                            p = p + 1
                        mass = string.atoi(self.smiles[pp:p+1])
            elif ch == "]":
                if not quoted:
                    # error, no opening ]
                    raise CEXsmilesError(smiles,p,"No opening ]")
                else:
                    quoted = 0
            elif ch == ".":
                self.fromat[lev] = None # disconnected parts
            # bond types
            elif ch == "=": bo = 2
            elif ch == "#": bo = 3
            elif ch == "-" and not quoted: bo = 1
            # atom charge
            elif ch == "-" or ch == "+":
                if not quoted:
                    # error charge not in []
                    raise CEXsmilesError(smiles,p,"Charge not in []")
                elif self.fromat[lev] is None:
                    # error charge precedes atomic symbol
                    raise CEXsmilesError(smiles,p,"Charge precedes atomic symbol")
                else:
                    charge = 0
                    sign = 1
                    if ch == "-": sign = -1
                    while self.smiles[p+1] in string.digits:
                        charge = 10*charge + string.atoi(self.smiles[p+1])
                        p = p + 1
                    if charge == 0: charge = 1
                    charge = sign*charge
                # allow for multiple + and - specifiers
                while self.smiles[p+1] == "+":
                    charge = charge + 1
                    p = p + 1
                while self.smiles[p+1] == "-":
                    charge = charge - 1
                    p = p + 1
                if charge != 0: self.SetFormalCharge(atom, charge)
            elif ch in string.digits or ch == "%" or ch == "^":
                # deal with ring closures
                if ch == "%":
                    if self.smiles[p+1] in string.digits and self.smiles[p+2] in string.digits:
                        ir = string.atoi(self.smiles[p+1:p+3])
                        p = p + 2
                    else:
                        # error expect 2 digits after %
                        raise CEXsmilesError(smiles,p,"Expect 2 digits after %")
                elif ch == "^":
                    if self.smiles[p+1] in string.digits and self.smiles[p+2] in string.digits and self.smiles[p+3] in string.digits:
                        ir = string.atoi(self.smiles[p+1:p+4])
                        p = p + 3
                    else:
                        #error expect 3 digits after ^
                        raise CEXsmilesError(smiles,p,"Expect 3 digits after ^")
                else:
                    ir = string.atoi(ch)
                if self.ringat[ir] is None:
                    self.ringat[ir] = self.fromat[lev]
                    self.ringbo[ir] = bo
                elif bo and self.ringbo[ir] and bo != self.ringbo[ir]:
                    #error conflicting closure bond orders
                    raise CEXsmilesError(smiles,p,"Conflicting closure bond orders")
                else:
                    if not bo: bo = 1
                    if self.ringbo[ir]: bo = self.ringbo[ir]
                    self.MakeBond(self.fromat[lev],self.ringat[ir],bo)
                    self.ringat[ir] = None
                    self.ringbo[ir] = 0
                bo = 0
            elif ch in "*ABCDEFGHIKLMNOPRSTUVWXYZ":
                # recognize atomic symbols
                atnum = -1
                if self.smiles[pp] in string.lowercase:
                    atnum = self.sym2num(self.smiles[p:p+2])
                if atnum > -1: p = p + 1
                else: atnum = self.sym2num(self.smiles[p])
                if atnum < 0:
                    #error bad atomic symbol
                    raise CEXsmilesError(smiles,p,"Bad atomic symbol")
                if not quoted and self.needquote(atnum):
                    # error symbol needs []'s
                    raise CEXsmilesError(smiles,p,"Symbol needs []")
                atom = self.MakeAtom(atnum)
                if not bo: bo = 1
                if (self.fromat[lev] is not None) and atom != self.fromat[lev]:
                    self.MakeBond(atom,self.fromat[lev],bo)
                self.fromat[lev] = atom
                if not quoted: imph = -1
                if mass > 0: self.SetAtomicMass(atom, mass)
                if quoted and atom is not None:
                    #deal with explict hydrogen counts
                    if self.smiles[p+1] != "H":
                        imph = 0
                    else:
                        imph = 1
                        p = p + 1
                        j = p
                        while self.smiles[p+1] in string.digits:
                            p = p + 1
                        if j < p: imph = string.atoi(self.smiles[j+1:p+1])
                if imph >= 0: self.SetHcount(atom,imph)
                # reset default attributes to undefined
                bo = 0
                charge = 0
                mass = 0
                imph = -1
            elif ch in string.whitespace:
                # extract molecul name from following text
                self.molname = self.smiles[p+1:-3]
                break
            elif ch == "\0":
                pass    #ignore guard characters
            else:
                # everything else is an error
                # error invalid character
                raise CEXsmilesError(smiles,p,"Invalid character")
            # end of while
            p = p + 1

class CEXprop:
    def __init__(self, tag, value):
        self.name = tag
        self.value = value
    def __str__(self):
        return self.name + "<" + self.value + ">"
        
class CEXchild(CEXprop):
    def __init__(self, tag, value):
        CEXprop.__init__(self, tag, value)
        self.proplist = []
    def __str__(self):
        str = self.name + "<" + self.value + ">"
        for p in self.properties():
            str = str + "\n" + p.__str__()
        return str
    def addProp(self, prop):
        self.proplist.append(prop)
    def properties(self):
        return self.proplist
class CEXroot(CEXchild):
    def __init__(self, tag, value):
        CEXchild.__init__(self, tag, value)
        self.childlist = []
    def addChild(self, child):
        self.childlist.append(child)
    def children(self):
        return self.childlist
    def __str__(self):
        str = self.name + "<" + self.value + ">"
        for p in self.properties():
            str = str + "\n" + p.__str__()
        for p in self.children():
            str = str + "\n" + p.__str__()
        return str
    
def readTree(cxstream):
    """Read tree of CEX object from stream"""
    (tag, value) = cxstream.readEntry()
    if not tag: return None
    if tag[0] != "$": return None
    root = CEXroot(tag,value)
    (tag, value) = cxstream.readEntry()
    if tag == None: return None
    while 1:
        if tag == "|": break
        if tag == None: break
        if tag[0] == "/":
            root.addProp(CEXprop(tag, value))
            (tag, value) = cxstream.readEntry()
        else:
            # Hardwired for root/child two level hierarchy
            child = CEXchild(tag, value)
            while 1:
                (tag, value) = cxstream.readEntry()
                if tag == "|": break
                if tag == None: break
                if tag[0] == "/":
                    child.addProp(CEXprop(tag, value))
                    continue
                else: break
            root.addChild(child)
    return root

def __follow_child(rec):
    print "  " + rec.name, rec.value
    for prop in rec.properties():
        print "    " + prop.name, prop.value

def spew(rec):
    print rec.name, rec.value
    for prop in rec.properties():
        print prop.name, prop.value
    for child in rec.children():
        __follow_child(child)

def selectChildren(rec, string):
    return filter(lambda x, string=string: x.name==string, rec.children())

def selectProperty(rec, string):
    for prop in rec.properties():
        if prop.name == string: return prop

if __name__ == "__main__":
    import StringIO
    def test(string):
        print "test: ",string
        s = StringIO.StringIO(string)
        c = CEXstream(s)
        print c.readEntry()
        s.close()
    test("|")
    test("tag<value>")
    test("  tag<value>")
    test("$tag<value>")
    test("/tag<value>")
    test("/tag_tag<value>")
    test('tag<"value">')
    test('tag<"value>">')
    test('tag<"""value>">')
    def test2(string):
        print "test2: ", string
        s = StringIO.StringIO(string)
        c = CEXstream(s)
        tree = readTree(c)
        spew(tree)
    test2("$root<test>|")
    test2("$root<test>/prop<value>|")
    test2("$root<test>child<value>|")
    test2("$root<test>/prop<value>/prop2<value2>|")
    test2("$root<test>/prop<value>/prop2<value2>child<valuec>|")
    test2("$root<test>/prop<value>/prop2<value2>child<valuec>/cprop<cv>|")
    def test2a(string):
        print "test2a: ", string
        s = StringIO.StringIO(string)
        c = CEXstream(s)
        tree = readTree(c)
        spew(tree)
        tree = readTree(c)
        spew(tree)
    test2a("$root<test>/prop<value>/prop2<value2>child<valuec>/cprop<cv>|$root2<test2>/prop<val>child<val>|")
    def test3(string):
        print "test3: ",string
        parser = CEXsmilesParser()
        try:
            parser.parse(string)
            print parser.molname
        except CEXsmilesError, data:
            print data

    test3("[C+2]")
    test3("[C++]")
    test3("[C+-]")
    test3("[C-2]")
    test3("[C--]")
    test3("[C-+]")
    test3("[CH3+2]")
    test3("N1#CC1")
    test3("N1#[CH3+2]C=1")
    test3("C%12CC%12")
    test3("C^123CC^123")
    test3("N1#[13CH3+2]C=1 test")
    test3("[N+1]C")
    test3("[N+]C")
    test3("N=[N+]=[N-]")
    test3("CC[[N]")
    test3("C=1CC-1")
    test3("[C]]")
    test3("C@1")
    test3("C+2")
    test3("[+2C]")
    test3("Si")
    test3("[Tx]")
    test3("C%1CC%1")
    test3("C^12CC^12")
    test3("[NH2+]")

