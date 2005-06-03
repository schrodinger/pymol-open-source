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

import chempy

def get_partial_charge(model):
    if chempy.feedback['verbose']:
        print ' '+str(__name__)+': generating partial charge keywords...'
    list = []
    c = -1
    for a in model.atom:
        list.append("CHARGE %d %6.4f\n" %(c,a.partial_charge))
        c = c - 1
    return list


def get_restrain_positions(model,flag,w_width,f_cnst):
    list = []
    n = 0
    c = 1
    mask = 1<<flag
    for a in model.atom:
        if (a.flags&mask):
            list.append("RESTRAIN-POSITION %5d %12.6f %12.6f %12.6f %6.3f %6.1f\n" %
                            (c,a.coord[0],a.coord[1],a.coord[2],w_width,f_cnst))
            n = n + 1
        c = c + 1
    if chempy.feedback['actions']:
        print ' '+str(__name__)+': %d atoms restrained using flag %d ...' % (n,flag)
        
    return list

def get_inactive(model,flag):
    list = []
    n = 0
    c = 1
    mask = 1<<flag
    for a in model.atom:
        if (a.flags&mask):
            list.append("INACTIVE %d\n" % (c))
            n = n + 1
        c = c + 1
    if chempy.feedback['actions']:
        print ' '+str(__name__)+': %d atoms fixed using flag %d ...' % (n,flag)
    return list

