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

try:
    import sglite

    from cmd import QuietException, \
          _feedback,fb_module,fb_mask 

    # xray.py 
    # This section contains python code for supporting
    # x-ray crystallography functions

    hex_to_rhom_xHM = {
        #  extended Hermann Mauguin symbol translation
        #  from CCP4's syminfo.lib
        
        'H 3'      : 'R 3: H',
        'H -3'     : 'R -3 :H',
        'H 3 2'    : 'R 3 2 :H',
        'H 3 m'    : 'R 3 m :H',
        'H 3 c'    : 'R 3 c :H',
        'H -3 2/m' : 'R -3 m :H',
        'H -3 m'   : 'R -3 m :H',
        'H -3 2/c' : 'R -3 c :H',
        'H -3 c'   : 'R -3 c :H',
        }
    
    def sg_sym_to_mat_list(sgsymbol):
        sgsymbol = hex_to_rhom_xHM.get(sgsymbol,sgsymbol)
        try:
            Symbols_Inp = sglite.SgSymbolLookup(sgsymbol)
            if(Symbols_Inp):
                HallSymbol = Symbols_Inp['Hall']
                SgOps = sglite.SgOps(HallSymbol)
                nLTr = SgOps.get_nLTr()
                fInv = SgOps.get_fInv()
                nSMx = SgOps.get_nSMx()
                result = []
                rb = float(sglite.SRBF)
                tb = float(sglite.STBF)
                for iLTr in xrange(nLTr):
                    for iInv in xrange(fInv):
                        for iSMx in xrange(nSMx):
                            Mx = SgOps.getLISMx(iLTr, iInv, iSMx, +1)
                            result.append([[ Mx[0]/rb, Mx[1]/rb, Mx[2]/rb, Mx[9 ]/tb],
                                           [ Mx[3]/rb, Mx[4]/rb, Mx[5]/rb, Mx[10]/tb],
                                           [ Mx[6]/rb, Mx[7]/rb, Mx[8]/rb, Mx[11]/tb],
                                           [        0,        0,        0,         1]] )
        except:
            if(_feedback(fb_module.symmetry,fb_mask.errors)):
                print "Symmetry-Error: Urecognized space group symbol '"+sgsymbol+"'."
            result = None
#        for a in result:
#            print 
#            for b in a:
#                print "%8.3f %8.3f %8.3f %8.3f"%(b[0],b[1],b[2],b[3])
        return result

except:
    print "Error: unable to import xray module (no symmetry support)."
    def sg_sym_to_mat_list(sgsymbol):
        return None
    
