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

from __future__ import print_function

import os
import shutil
import glob
import sys

from chempy import feedback

# "do" is the preferred command for running tinker

def do(command,in_prefix,run_prefix,out_prefix,tokens,capture=None):
    if feedback['tinker']:
        print(" "+str(__name__)+': creating temporary files "%s.*"' % (run_prefix))
        print(" "+str(__name__)+': launching %s...' % command)
        c = 1
        for a in tokens:
            print(" "+str(__name__)+': input %d = %s' % (c,a))
            c = c + 1
        if hasattr(sys.stdout,"flush"):
            sys.stdout.flush()
    for a in glob.glob(run_prefix+".*"):
        os.unlink(a)
    for a in glob.glob(out_prefix+".*"):
        os.unlink(a)
    for src in glob.glob(in_prefix+".*"):
        dst = src.split('.')
        dst = run_prefix+'.'+dst[len(dst)-1]
        shutil.copyfile(src,dst)
    if capture==1:
        pipe = os.popen(bin_path+command+"> "+run_prefix+".out","w")
    elif capture==2:
        pipe = os.popen(bin_path+command+" | tee "+run_prefix+".out","w")
    else:
        pipe = os.popen(bin_path+command,"w")
    if not pipe:
        print("Error: can't run tinker!!!")
        raise RunError
    for a in tokens:
        pipe.write(a+"\n")
    pipe.close()
# NFS workaround (flushes the directory cache so that glob will work)
    try: os.unlink(".sync")
    except: pass
    f = open(".sync",'w')
    f.close()
#
    for src in glob.glob(run_prefix+".*_2"):
        dst = src.replace('_2','')
        if os.path.exists(dst):
            os.unlink(dst)
#      os.rename(src,dst)    rename can fail over NFS (remote action)
        shutil.copyfile(src,dst)
# sloppy workaround for buggy NFS on linux
        os.unlink(src)
    for src in glob.glob(run_prefix+".*"):
        dst = src.split('.')
        dst = out_prefix+'.'+dst[len(dst)-1]
        if os.path.exists(dst):
            os.unlink(dst)
#      os.rename(src,dst)    rename can fail over NFS (remote action)
        shutil.copy(src,dst)
        os.unlink(src)
    for a in glob.glob(in_prefix+".*"):
        os.unlink(a)
    if feedback['tinker']:
        print(" "+str(__name__)+': %s job complete. ' % command)
        print(" "+str(__name__)+': creating output files "%s.*"' % (out_prefix))

#  DEPRECATED

prefix = "tinker_run"

def run(command,in_prefix,out_prefix,tokens,capture=None):
    if feedback['tinker']:
        print(" "+str(__name__)+': creating temporary files "%s.*"' % (prefix))
        print(" "+str(__name__)+': launching %s...' % command)
        c = 1
        for a in tokens:
            print(" "+str(__name__)+': input %d = %s' % (c,a))
            c = c + 1
        if hasattr(sys.stdout,"flush"):
            sys.stdout.flush()
    for a in glob.glob(prefix+".*"):
        os.unlink(a)
    for a in glob.glob(out_prefix+".*"):
        os.unlink(a)
    for src in glob.glob(in_prefix+".*"):
        dst = src.split('.')
        dst = prefix+'.'+dst[len(dst)-1]
        shutil.copyfile(src,dst)
    if capture:
        pipe = os.popen(bin_path+command+"> "+out_prefix+".out","w")
    else:
        pipe = os.popen(bin_path+command,"w")
    if not pipe:
        print("Error: can't run tinker!!!")
        raise RunError
    for a in tokens:
        pipe.write(a+"\n")
    pipe.close()
    for src in glob.glob(prefix+".*_2"):
        dst = src.replace('_2','')
        if os.path.exists(dst):
            os.unlink(dst)
#      os.rename(src,dst)    rename can fail over NFS (remote action)
        shutil.copy(src,dst)
        os.unlink(src)
    for src in glob.glob(prefix+".*"):
        dst = src.split('.')
        dst = out_prefix+'.'+dst[len(dst)-1]
        if os.path.exists(dst):
            os.unlink(dst)
#      os.rename(src,dst)    rename can fail over NFS (remote action)
        shutil.copy(src,dst)
    if feedback['tinker']:
        print(" "+str(__name__)+': %s job complete. ' % command)
        print(" "+str(__name__)+': creating output files "%s.*"' % (out_prefix))

if 'TINKER_PATH' in os.environ:
    base = os.environ['TINKER_PATH']
    bin_path = base + '/bin/'
    params_path = base + '/params/'
elif 'FREEMOL_ETC' in os.environ:
    base = os.environ['FREEMOL_ETC'] + '/tinker'
    bin_path = base + '/bin/'
    params_path = base + '/params/'
else:
    base = ''
    bin_path = ''
    params_path = ''

if 'PYMOL_PATH' in os.environ:
    pymol_path = os.environ['PYMOL_PATH']
    test_path = pymol_path + '/data/chempy/tinker/'
    if os.path.exists(test_path):
        params_path = test_path
