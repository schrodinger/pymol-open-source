
dali_file = "dali.txt"

pdb_dir = "dali_pdb"

max_pairs = 10

def fetch_pdb(pdbCode,outFile):
    import urllib
    import gzip
    import os
    import string

    remoteCode = string.upper(pdbCode)
    if not os.path.exists(pdb_dir):
        os.mkdir(pdb_dir)
    if not os.path.exists(outFile):
        try:
            filename = urllib.urlretrieve(
                'http://www.rcsb.org/pdb/cgi/export.cgi/' +
                remoteCode + '.pdb.gz?format=PDB&pdbId=' +
                remoteCode + '&compression=gz')[0]
        except:
            print "warning: %s not found.\n"%pdbCode
        else:
            if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
                try:
                    abort = 0
                    open(outFile, 'w').write(gzip.open(filename).read())
                    print "fetched: %s"%(pdbCode)
                except IOError:
                    abort = 1
                if abort:
                    os.remove(outFile)
            else:
                print "warning: %s not valid.\n"%pdbCode
            os.remove(filename)

from pymol import cmd
from string import strip
import os

seen = {}
input = open(dali_file).readlines()
input_state = 0
while 1:
    try:
        line = input.pop(0)
    except IndexError:
        break
    if input_state == 0:
        if line[0:11]=='## MATRICES':
            line = input.pop(0)
            if line[0:12]=='  NR. STRID1':
                input_state = 1
    elif input_state == 1:
        if strip(line)=='':
            input_state = 2
        elif line[4:5]==':':
            trg = strip(line[6:12])
            src = strip(line[13:19])
            trg_code = trg[0:4]
            src_code = src[0:4]
            if not seen.has_key(trg_code):
                trg_file = pdb_dir+os.sep+trg_code+".pdb"
                fetch_pdb(trg_code,trg_file)
                cmd.load(trg_file)
                seen[trg_code]=1
            if not seen.has_key(src_code):
                src_file = pdb_dir+os.sep+src_code+".pdb"
                fetch_pdb(src_code,src_file)
                cmd.load(src_file)
                seen[src_code]=1
            matrix = []
            for a in range(0,3):
                matrix.append(float(strip(line[29:38])))
                matrix.append(float(strip(line[39:48])))
                matrix.append(float(strip(line[49:58])))
                matrix.append(float(strip(line[59:78])))
            matrix.extend([0.0,0.0,0.0,1.0])
            max_pairs = max_pairs - 1
    if max_pairs<0:
        break
            
