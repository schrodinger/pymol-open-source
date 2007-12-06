#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

import string
import re
import copy

from chempy import io

# OMG what a horrid file format.

single_quote_re = re.compile(r"'[^']*[']*'") # doesn't yet handle ESC
double_quote_re = re.compile(r'"[^"]*["]*"') # ditto
bracket_quote_re = re.compile(r'\[[^\]]*\]') # ditto

class CIFRec:

    def get_quoted_value(self):
        if self.line[0:1] == "'":
            mo = single_quote_re.match(self.line)
        elif self.line[0:1] == '"':
            mo = double_quote_re.match(self.line)
        elif self.line[0:1] == '[':
            mo = bracket_quote_re.match(self.line)
        else:
            mo = None
        if mo != None:
            result = self.line[:mo.end()]
            self.line = self.line[mo.end():]
            if len(string.strip(self.line[0:1])): # followed by non-whitespace...
                self.line = result[-1:] + self.line
                result = result[:-1] + self.get_quoted_value()
            self.check_line()
        else:
            result = None # shouldn't happen
        return result
    
    def get_delimited_value(self):
        result = ''
        while 1:
            if self.line[0:1] == ';':
                self.line = self.line[1:]
                self.check_line()
                break
            result = result + self.line
            self.next_line()
        return result

    def get_next_word(self):
        result = None
        mo = re.search("\s+", self.line)
        if mo == None:
            result = self.line
            self.next_line()
        else:
            result = self.line[:mo.start()]
            self.line = self.line[mo.end():]
            self.check_line()
        return result
    
    def trim_leading_whitespace(self):
        while self.line != None:
            mo = re.match("\s",self.line)
            if mo != None:
                self.line = self.line[mo.end():]
                self.check_line()
            else:
                break

    def get_next_value(self):
        self.trim_leading_whitespace()
        if self.line != None:
            if self.line[0:1] == '_':
                return None
            elif self.line[0:1] == ';':
                self.line = self.line[1:]
                self.check_line()
                return self.get_delimited_value()
            elif string.lower(self.line[0:5]) in ['loop_','data_','save_']:
                return None
            elif string.lower(self.line[0:6]) == 'GLOBAL':
                return None
            elif self.line[0:1] in [ "'", '"']: #  '[' ]: bracket quote fubar with PDB's mmCIF data???
                return self.get_quoted_value()
            else:
                return self.get_next_word()
        return None
        
    def check_line(self):
        if self.line != None:
            while not len(string.strip(self.line)):
                self.next_line()
                if self.line == None:
                    break

    def next_line(self):
        self.line = None
        while len(self.list):
            self.line = self.list.pop()
            # nuke hash comments 
            hash = string.find(self.line,'#')
            if hash>=0:
                self.line = self.line[0:hash]
            # and only return non-blank lines
            if len(string.strip(self.line)):
                break
#        print "next_line:", self.line,
        
    def parse_loop_body(self,fields):
        fields_rev = fields.reverse()
        len_fields = len(fields)
        records = []
        record = []
        cnt = len_fields
        while self.line != None: 
            value = self.get_next_value()
            if value == None:
                break
            else:
                cnt = cnt - 1
                if not cnt:
                    cnt = len_fields
                    if len(record) == len_fields:
                        records.append(record)
                    record = []
            record.append(value)
#                print "loop_read [%s]=[%s]"%(fields[0],value)
        if len(record) == len_fields:
            records.append(record)
        self.loops.append( (fields,records) )
        
    def parse_loop(self):
#        print "parsing loop..."
        fields = []
        while self.line != None: 
            if string.lower(self.line[0:5])=='loop_':
                break
            elif self.line[0:1]=='_':
                fields.append(string.lower(string.strip(self.line)))
            else:
                self.parse_loop_body(fields)
                break
            self.next_line()
                
    def parse_name_value(self):
        data_name = self.get_next_word()
        data_value = self.get_next_value()
        self.key_value[data_name] = data_value
#        print "data_read [%s]=[%s]"%(data_name,data_value)
        
    def parse_normal(self):
        while self.line != None:
            if self.line[0:1] == '_': # data name
                self.parse_name_value()
            elif string.lower(self.line[0:5])=='loop_':
                self.line = self.line[5:]
                self.check_line()
                self.parse_loop()
            else: # shouldn't happen
                print "unhandled: [%s]"%self.line
                self.next_line()
            
    def __init__(self,cif_list):
        cif_list.reverse()
        self.list = cif_list
        data_line = self.list.pop()
        self.loops = []
        self.key_value = {}
        self.data_name = string.strip(data_line[5:])
        self.next_line()
        self.parse_normal()
        print ' CIF: For data block "%s"...'%self.data_name
        print " CIF: Read %d key/value pair(s)."%len(self.key_value)
        print " CIF: Read %d table(s)."%len(self.loops)
        # does this record have 
    def toList(self):
        return r

class CIF:
    
    def __init__(*args):
        mode = 'r'
        if len(args)<2:
            raise ValueError
        self = args[0]
        self.input_line = None
        fname = args[1]
        if len(args)==3:
            mode = args[2]
        self.mode = mode
        self.at_eof = 0
        if mode not in ('w','r','wa','pf','url'):
            print " CIF: bad mode"
            return None
        if mode=='pf': # pseudofile
            self.file = fname
        elif (mode[0:1]=='r') and (string.find(fname,':')>1):
            # does this look like a URL? (but not a DOS path)
            from urllib import urlopen
            self.file = urlopen(fname)
        else:
            self.file = open(fname,mode)

    def write(self,rec):
        lst = rec.toList()
        for a in lst:
            self.file.write(a)
        self.file.write('$$$$\n')
        
    def read(self): # returns CIFRec or None at end of file
        cur = []
        while 1:
            if self.input_line == None:
                s = self.file.readline()
            else:
                s = self.input_line
                self.input_line = None
            if not s: # end of file
                if len(cur)>0:
                    return CIFRec(cur)
                else:
                    return None
            elif string.lower(s[0:5]) == r'data_': # signals start of new record
                if len(cur)>0:
                    self.input_line = s # save for next time
                    return CIFRec(cur)
                else:
                    cur.append(s)
            else:
                cur.append(s)
    def close(self):
        self.file.close()
        
    
