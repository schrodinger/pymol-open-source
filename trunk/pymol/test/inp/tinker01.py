# python

from chempy import io
from chempy.tinker import now
from chempy import tinker
from chempy.tinker import keyword
from chempy import protein
import os

a = io.pdb.fromFile("dat/il2.pdb")

a=protein.generate(a)

inp_prefix = 'tinker_inp'
out_prefix = 'tinker_out'
def_params = 'chempy.prm'

kw = [
   "parameters        "+tinker.params_path+def_params+"\n",
   "lights\n",
   "verbose\n",
   "randomseed         1234567890\n",
   "cutoff             6.0\n", 
   "overwrite\n",
   ]

io.xyz.toFile(a,inp_prefix+".xyz")

kw.extend(keyword.get_partial_charge(a))   

f=open(inp_prefix+".key",'w')
f.writelines(kw)
f.close()

tinker.run('dynamic',inp_prefix,out_prefix,
              [tinker.prefix,
					'10',
               '1.0',
					'0.01',
               '300'],capture=1)

io.xyz.updateFromFile(a,out_prefix+".001")

c = 0
for at in a.atom:
   if c > 500: break
   print "%-4s %-4s %12.6f %12.6f %12.6f" % (
      at.name,at.resi,at.coord[0],at.coord[1],at.coord[2])
   c = c + 1
   
os.system("touch .no_fail tinker_*")
os.system("/bin/rm .no_fail tinker_*")

