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

import copy

def __list_sum(the_list):
    return reduce(lambda x,y:x+y,the_list)

def combine_fragments(*arg,**kw):
    '''
    merge_fragments(target,source1,source2,source3,...)
    
    WARNING: This is not how you *should* combine fragments,
    rather this is just a convenient kludge for getting
    something remotely realistic in a big hurry.
    
    NOTE: atom names must be unique throughout -- especially
    including the "capping" atoms used during the ab-initio
    calculations (you can use PyMOL to rename these if they
    happen to coincide).
'''
    
    if kw.has_key('net_charge'):
        net_charge = kw['net_charge']
    else:
        net_charge = None
    if len(arg)<1:
        raise TypeError('invalid arguments')
    dst = copy.deepcopy(arg[0])
    frg_lst = arg[1:]
    n_frg = len(frg_lst)
    n_dst_atm = len(dst.atom)

    # create dictionary/indices for shared atoms
    
    n_tot = 0
    dst_dict = {}
    for a in dst.atom:
        sig = a.name
        dst_dict[sig] = n_tot
        n_tot = n_tot + 1

    # create fragment membership list for destination atoms
    
    members = []
    for a in range(n_tot):
        members.append([])
    
    # create indices for unique atoms

    for fragment in frg_lst:
        c = 0
        for a in fragment.atom:
            sig = a.name
            if dst_dict.has_key(sig):
                a.chg_index = dst_dict[sig]
                members[a.chg_index].append(fragment)
            else:
                # find an attached atom that *is* in the destination structure
                attached_index = None
                for b in fragment.bond:
                    attached = None
                    if b.index[0]==c:
                        attached = fragment.atom[b.index[1]]
                    elif b.index[1]==c:
                        attached = fragment.atom[b.index[0]]
                    if attached!=None:
                        if dst_dict.has_key(attached.name):
                            attached_index = dst_dict[attached.name]
                            break
                a.chg_index = n_tot
                if attached_index!=None:
                    a.attached_index = attached_index
                n_tot = n_tot + 1
            c = c + 1
                
    # create array for charges
    
    chg = []
    for a in range(n_frg):
        chg.append([0.0] * n_tot)
    
    # now load and count measurements

    cnt = [0] * n_tot
    c = 0
    for fragment in frg_lst:
        for a in fragment.atom:
            index = a.chg_index
            if (index < n_dst_atm):
                chg[c][index] = chg[c][index] + a.partial_charge
                cnt[index] = cnt[index] + 1
            else:
                aindex = a.attached_index
                chg[c][aindex] = chg[c][aindex] + a.partial_charge
        c = c + 1

    # now compute average charges for each destination atom

    avg = []
    for index in range(n_dst_atm):
        if cnt[index]:
            tmp_lst = []
            for a in chg:
                tmp_lst.append(a[index])
            avg.append(__list_sum(tmp_lst)/cnt[index])
        else:
            avg.append(0.0)
        
    # correct total charge

    chg_sum = __list_sum(avg)

    print "chg_sum",chg_sum
    if net_charge != None:
        chg_diff = net_charge - chg_sum
        chg_adjust = chg_diff / n_dst_atm
    else:
        chg_adjust = 0.0
        
    c = 0
    for a in dst.atom:
        a.partial_charge = avg[c] + chg_adjust
        c = c + 1

    return dst

