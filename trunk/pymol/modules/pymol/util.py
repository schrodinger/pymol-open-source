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

import cmd
import math
import string
import glob
import pymol

def cbag(selection):
   s = str(selection)
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("carbon","(elem C and "+s+")")

def cbac(selection):
   s = str(selection)
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("cyan","(elem C and "+s+")")

def cbay(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("yellow","(elem C and "+s+")")

def cbas(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("salmon","(elem C and "+s+")")

def cbap(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("purple","(elem C and "+s+")")

def cbaw(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("hydrogen","(elem C and "+s+")")

def cbab(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("slate","(elem C and "+s+")")

def mrock(first,last,angle=30,phase=0,loop=1,axis='y'):
   first=int(first)
   last=int(last)
   angle=float(angle)
   phase=float(phase)
   loop=int(loop)
   nstep = (last-first)+1
   if nstep<0:
      nstep = 1
   if loop:
      subdiv = nstep
   else:
      subdiv = nstep+1
   ang_cur = math.pi*phase/180
   ang_inc = 2*math.pi/subdiv
   ang_cur = ang_cur - ang_inc
   a = 0
   while a<nstep:
      last = angle*math.sin(ang_cur)
      ang_cur = ang_cur + ang_inc
      disp = angle*math.sin(ang_cur)
      diff = disp-last
      com = "mdo %d:turn %s,%8.3f" % (first+a,axis,diff)
      cmd.do(com)
      a = a + 1

def mroll(first,last,loop=1,axis='y'):
   first=int(first)
   last=int(last)
   loop=int(loop)
   n = last - first
   if loop:
      step = 2*math.pi/(n+1)
   else:
      step = 2*math.pi/n   
   a = 0
   deg = (180*step/math.pi)
   while a<=n:
      com = "mdo %d:turn %s,%8.3f" % (first+a,axis,deg)
      cmd.do(com)
      a = a + 1

def hbond(a,b,cutoff=3.3):
   st = "(%s and (%s around %4.2f) and elem N,O),(%s and (%s around %4.2f) and elem N,O),%4.2f" % (a,b,cutoff,b,a,cutoff,cutoff)
   cmd.dist("hbond",st)
        
def mload(*args):
   nam = "mov"
   if len(args)>1:
      nam = args[1]
   fils = glob.glob(args[0])
   fils.sort()
   if not len(fils):
      print "Error: no matching files"
   else:
      for a in fils:
         cmd.load(a,nam)
   
def cbc(selection='(all)'): # NOT THREAD SAFE
   '''
   Color all chains a different color
   '''
   pymol.stored.chain = {}
   cmd.iterate("(%s)"%selection,"stored.chain[chain]=1")
   c = 7
   for a in pymol.stored.chain.keys():
      if len(a):
         print ("%d,(chain %s)"%(c,a))
         cmd.color("%d"%c,"(chain %s)"%a)
         c = c + 1

color_chains = cbc

def sum_charge(*arg): # NOT THREAD SAFE
   result = None
   try:
      obj = "all"
      if len(arg):
         obj = arg

      pymol.stored._sum_charge = 0.0
      cmd.iterate("(%s)"%obj,
                  "stored._sum_charge=stored._sum_charge+partial_charge")
      result = pymol.stored._sum_charge
      print " sum_charge: %6.4f"%result
   except:
      print " sum_charge: an error occurred."
   return result


def ff_copy(src,dst): # NOT THREAD SAFE
   pymol._rcopy = pymol.Scratch_Storage()
   pymol._rcopy.pc={}
   pymol._rcopy.tt={}
   cmd.iterate("(%s)"%src,"_rcopy.pc[name]=partial_charge")
   cmd.alter("(%s)"%dst,"partial_charge=_rcopy.pc[name]")
   cmd.iterate("(%s)"%src,"_rcopy.tt[name]=text_type")
   cmd.alter("(%s)"%dst,"text_type=_rcopy.tt[name]")
   del pymol._rcopy
   
def b2vdw(*arg):
   if not len(arg):
      sele = 'all'
   else:
      sele = arg[0]
   # use B values to create RMS VDW spheres
   # rms = sqrt(b/(8*(PI^2)))
   cmd.alter("(%s)"%sele,"vdw=math.sqrt(b/78.9568352087)")
   
def phipsi(selection="(pk1)"): # NOT THREAD SAFE
   n_sele =   "((byres (%s)) & name n)"%selection
   c_sele =   "((byres (%s)) & name c)"%selection
   ca_sele =  "((byres (%s)) & name ca)"%selection
   cm_sele = "((neighbor (%s)) and not (byres (%s)))"%(n_sele,n_sele)
   np_sele = "((neighbor (%s)) and not (byres (%s)))"%(c_sele,c_sele)
   cmd.feedback("push")
   cmd.feedback("disable","selector","everythin")
   cm_cnt = cmd.select("pp_cm",cm_sele)
   n_cnt = cmd.select("pp_n",n_sele)
   c_cnt = cmd.select("pp_c",c_sele)
   ca_cnt = cmd.select("pp_ca",ca_sele)
   np_cnt = cmd.select("pp_np",np_sele)
   if(cm_cnt and n_cnt and ca_cnt and c_cnt):
      phi = cmd.get_dihedral("pp_c","pp_ca","pp_n","pp_cm")
   else:
      phi = None
   if(n_cnt and ca_cnt and c_cnt and np_cnt):
      psi = cmd.get_dihedral("pp_np","pp_c","pp_ca","pp_n")
   else:
      psi = None
   cmd.feedback("pop")
   return (phi,psi)

def rainbow(selection="(name ca and alt '',A)"): # NOT THREAD SAFE

   cmd.feedback("push")
   cmd.feedback("disable","executive","actions")

   # your basic rainbow...
   
   list = [(255,0,0),(255,128,0),
           (255,255,0),(0,255,0),
           (0,255,255),(0,128,255),
           (0,0,255)]
   #
   last = list.pop(0)
   cmd.set_color("rainbow000",[last[0]/255.0,last[1]/255.0,last[2]/255.0])
   c = 1
   for a in list:
      for b in range(1,21):
         b0 = b/20.0
         b1 = 1.0-b0
         cname = "rainbow%03d"%c
         r = last[0]*b1+a[0]*b0
         g = last[1]*b1+a[1]*b0
         b = last[2]*b1+a[2]*b0
         cmd.set_color(cname,[r/255.0,g/255.0,b/255.0])
         c = c + 1
      last = a

   cas = cmd.index("((byres ("+selection+")) and name ca and not het)")
   l = len(cas)
   if not len(cas):
      return
   c = 0
   for a in cas:
      col = int((120*c)/l)
      cmd.color("rainbow%03d"%col,"(byres %s`%d)"%a)
      c = c + 1

   cmd.feedback("pop")
   
def ss(selection="(name ca and alt '',A)"): # NOT THREAD SAFE

   cmd.feedback("push")
   cmd.feedback("disable","executive","actions")
   
   ss_pref = "sss"
   sss1 = ss_pref+"1"
   cnt = cmd.select(sss1,"((byres ("+selection+")) and name ca and not het)")
   print " util.ss: initiating secondary structure assignment on %d residues."%cnt
   cas = cmd.index(sss1)
   if not len(cas):
      return

   cmd.cartoon("auto",sss1)

   print " util.ss: extracting residue sequence..."

   # get CA list
   
   res_list = []
   pymol._ss = pymol.Scratch_Storage()
   pymol._ss.res_list = res_list
   cmd.iterate(sss1,'_ss.res_list.append((model,index))')

   # generate atom-to-residue conversion dictionaries

   n_dict = {}
   o_dict = {}
   scr_dict = {}
   pymol._ss.n_dict = n_dict
   pymol._ss.o_dict = o_dict
   pymol._ss.scr_dict = scr_dict
   cmd.iterate(sss1,'_ss.scr_dict[(model,index)]=(segi,chain,resi)')
   cmd.iterate("((byres "+sss1+") and n;n)",'_ss.n_dict[(segi,chain,resi)] = (model,index)')
   cmd.iterate("((byres "+sss1+") and n;o)",'_ss.o_dict[(segi,chain,resi)] = (model,index)')

   # find gaps

   gap_list = [None,None,None,None]  # large enough to distinguish i+4 interations from gaps
   last = None
   for a in res_list:
      if last!=None:
         if(cmd.count_atoms(
            "((neighbor(neighbor(neighbor (%s`%d)))) and (%s`%d))"%
            (last[0],last[1],a[0],a[1]),quiet=1)==0):
            gap_list.extend([None,None,None,None])
         gap_list.append(a)
      last = a
   gap_list.extend([None,None,None,None])

   print " util.ss: analyzing phi/psi angles..."

   ss = {}
   ss[None]=None
   
   # intialize all residues to loop
   
   for a in cas:
      ss[a] = 'L'
      
   # make decisions based on phi/psi

   for a in cas:
      (phi,psi) = phipsi("(%s`%d)"%(a[0],a[1]))
      if (phi!=None) and (psi!=None):
         if ((phi<-45) and (phi>-160) and (psi>90)): # beta?
            ss[a]='S'
         elif ((phi<-45) and (phi>-160) and
               (psi>-80) and (psi<-30)): # helix?
            ss[a]='H'

   # find all pairwise hydrogen bonds and make note of them in dict

   hb = cmd.find_pairs("((byres "+sss1+") and n;n)",
                       "((byres "+sss1+") and n;o)",mode=1,angle=35)
   hb_dict = {}
   for a in hb:
      hb_dict[a] = 1
      print a

   # check to insure that all helical residues have at least an i +/- 4
   # hydrogen bond

   for c in xrange(4,len(cas)-4):
      a = cas[c]
      if ss[a]=='H':
         aN = n_dict[scr_dict[a]]
         aO = o_dict[scr_dict[a]]
         am4O = o_dict[scr_dict[cas[c-4]]]
         ap4N = n_dict[scr_dict[cas[c+4]]]
         if not hb_dict.has_key((aN,am4O)):
            if not hb_dict.has_key((ap4N,aO)):
               ss[a]='L'

   # extend helices if hydrogen bonding requirements are met

   repeat = 1
   while repeat:
      repeat = 0
      for c in xrange(4,len(cas)-4):
         if ss[cas[c+1]]=='H': 
            a = cas[c]
            if ss[a]!='H': # N-terminal end
               aO = o_dict[scr_dict[a]]
               ap4N = n_dict[scr_dict[cas[c+4]]]
               if hb_dict.has_key((ap4N,aO)):
                  ss[a]='H'
                  repeat = 1
         if ss[cas[c-1]]=='H':
            a = cas[c]
            if ss[a]!='H': # C-terminal end
               aN = n_dict[scr_dict[a]]
               am4O = o_dict[scr_dict[cas[c-4]]]
               if hb_dict.has_key((aN,am4O)):
                  ss[a]='H'
                  repeat = 1

   # check to insure that all beta residues have at least a pair
   # of interacting pairs

   if 0:
      for c in xrange(4,len(cas)-4):
         a = cas[c]
         if ss[a]=='S':
            print c,cas[c]
            aN = n_dict[scr1]
            aO = o_dict[scr1]
            aS = scr_dict[a]
            aSp1 = cas[c+1]
            aSm1 = cas[c-1]
            found = 0
            # antiparallel, Na->Ob, Nb+2->Oa-2
            if not found:

               am1O = o_dict[scr_dict[cas[c-4]]]
               ap1N = n_dict[scr_dict[cas[c+4]]]
               if not hb_dict.has_key((aN,am4O)):
                  if not hb_dict.has_key((ap4N,aO)):
                     ss[a]='L'

   
   # assign protein
   
   cmd.alter(sss1,"ss ='L'")
   for a in cas:
      if ss[a]!='L':
         cmd.alter("(%s`%d)"%a,"ss='%s'"%ss[a])

   cmd.feedback("pop")

   del pymol._ss # IMPORTANT
   
   #
#   print conn_hash.keys()
   print " util.ss: assignment complete."


