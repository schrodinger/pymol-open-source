
from chempy import tinker
from chempy import io
from chempy.tinker import keyword

inp_prefix = 'tinker_inp'
out_prefix = 'tinker_out'
def_params = 'chempy.prm'

def_keywords = [
   "parameters        "+tinker.params_path+def_params+"\n",
   "lights\n",
   "verbose\n",
   "randomseed         1234567890\n",
   "cutoff             8.0\n",  # very small
   ]

def _write_keys(list):
   f=open(inp_prefix+".key",'w')
   f.writelines(list)
   f.close()
   
def minimize(model,gradient=5.0,maxiter=None,kw=None):
   io.xyz.toFile(model,inp_prefix+".xyz")
   if not kw:
      kw = []
      kw.extend(def_keywords)
   if maxiter:
      kw.append("maxiter %d\n" % maxiter)
   kw.append("overwrite\n")
   kw.extend(keyword.get_partial_charge(model))   
   _write_keys(kw)
   tinker.run('minimize',inp_prefix,out_prefix,[tinker.prefix,str(gradient)])
   io.xyz.updateFromFile(model,out_prefix+".xyz")

def dynamics(model,steps=100,timestep=1.0,dumps=0.1,temperature=300,kw=None):
   io.xyz.toFile(model,inp_prefix+".xyz")
   if not kw:
      kw = []
      kw.extend(def_keywords)
   kw.append("overwrite\n")
   kw.extend(keyword.get_partial_charge(model))   
   _write_keys(kw)
   tinker.run('dynamic',inp_prefix,out_prefix,
              [tinker.prefix,str(steps),
               str(timestep),str(dumps),
               str(temperature)])
   io.xyz.updateFromFile(model,out_prefix+".001")


def analyze(model,gradient=5.0,maxiter=None,kw=None):
   io.xyz.toFile(model,inp_prefix+".xyz")
   if not kw:
      kw = []
      kw.extend(def_keywords)
   if maxiter:
      kw.append("maxiter %d\n" % maxiter)
   kw.append("overwrite\n")
   kw.extend(keyword.get_partial_charge(model))   
   _write_keys(kw)
   tinker.run('analyze',inp_prefix,out_prefix,[tinker.prefix,"E"])
   io.xyz.updateFromFile(model,out_prefix+".xyz")


