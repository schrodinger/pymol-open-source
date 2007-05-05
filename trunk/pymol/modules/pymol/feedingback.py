
import sys
import cmd
from cmd import Shortcut, is_string
from cmd import fb_module, fb_mask, fb_action

import _cmd
import string

def _feedback(module,mask,_self=cmd): # feedback query routine
    # WARNING: internal routine, subject to change      
    r = 0
    module = int(module)
    mask = int(mask)
    if module>0:
        try:
            _self.lock(_self)
            r = _cmd.feedback(_self._COb,module,mask)
        finally:
            _self.unlock(-1)
    else:
        if fb_dict.has_key(module):
            r = fb_dict[module]&mask
    return r


fb_action_sc = Shortcut(fb_action.__dict__.keys())

fb_module_sc = Shortcut(fb_module.__dict__.keys())


fb_mask_sc = Shortcut(fb_mask.__dict__.keys())

fb_dict ={}

for a in fb_module.__dict__.keys():
    vl = getattr(fb_module,a)
    if vl<0:
        fb_dict[vl] = 0x1F # default mask

fb_debug = sys.stderr # can redirect python debugging output elsewhere if desred...

def feedback(action="?",module="?",mask="?",_self=cmd):
    '''
DESCRIPTION

    "feedback" allows you to change the amount of information output by pymol.

USAGE

    feedback action,module,mask

    action is one of ['set','enable','disable']
    module is a space-separated list of strings or simply "all"
    mask is a space-separated list of strings or simply "everything"

NOTES:

    "feedback" alone will print a list of the available module choices

PYMOL API

    cmd.feedback(string action,string module,string mask)

EXAMPLES

    feedback enable, all , debugging
    feedback disable, selector, warnings actions
    feedback enable, main, blather

DEVELOPMENT TO DO

    Add a way of querying the current feedback settings.
    Check C source code to make source correct modules are being used.
    Check C source code to insure that all output is properly
    Update Python API and C source code to use "quiet" parameter as well.
        '''
    r = None

    # validate action

    if action=="?":
        print " feedback: possible actions: \nset, enable, disable"
        act_int = 0
    else:
        act_kee = fb_action_sc.interpret(action)
        if act_kee == None:
            print "Error: invalid feedback action '%s'."%action
            if _raising():
                raise QuietException
            else:
                return None
        elif not is_string(act_kee):
            print "Error: ambiguous feedback action '%s'."%action
            print action_amb
            if _raising():
                raise QuietException
            else:
                return None
        act_int = int(getattr(fb_action,act_kee))

    if (act_int<3) and ("?" in [action,module,mask]):
        if module=="?":
            print " feedback: Please specify module names:"
            lst = fb_module.__dict__.keys()
            lst.sort()
            for a in lst:
                if a[0]!='_':
                    print "   ",a
        if mask=="?":
            print " feedback: Please specify masks:"
            lst = fb_mask.__dict__.keys()
            lst.sort()
            for a in lst:
                if a[0]!='_':
                    print "   ",a
    else:
        if (act_int>=3):
            module='all'
            mask='everything'

        # validate and combine masks

        mask_int = 0
        mask_lst = string.split(mask)
        for mask in mask_lst:
            mask_kee = fb_mask_sc.interpret(mask)
            if mask_kee == None:
                print "Error: invalid feedback mask '%s'."%mask
                if _raising(): raise QuietException
                else: return None
            elif not is_string(mask_kee):
                print "Error: ambiguous feedback mask '%s'."%mask
                if _raising(): raise QuietException
                else: return None
            mask_int = int(getattr(fb_mask,mask_kee))

        # validate and iterate modules

        mod_lst = string.split(module)
        for module in mod_lst:
            mod_kee = fb_module_sc.interpret(module)
            if mod_kee == None:
                print "Error: invalid feedback module '%s'."%module
                if _raising(): raise QuietException
                else: return None
            elif not is_string(mod_kee):
                print "Error: ambiguous feedback module '%s'."%module
                if _raising(): raise QuietException
                else: return None
            mod_int = int(getattr(fb_module,mod_kee))
            if mod_int>=0:
                try:
                    _self.lock(_self)
                    r = _cmd.set_feedback(_self._COb,act_int,mod_int,mask_int)
                finally:
                    _self.unlock(_self=_self)
            if mod_int<=0:
                if mod_int:
                    if act_int==0:
                        fb_dict[mod_int] = mask_int
                    elif act_int==1:
                        fb_dict[mod_int] = fb_dict[mod_int] | mask_int
                    elif act_int==2:
                        fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
                else:
                    for mod_int in fb_dict.keys():
                        if act_int==0:
                            fb_dict[mod_int] = mask_int
                        elif act_int==1:
                            fb_dict[mod_int] = fb_dict[mod_int] | mask_int
                        elif act_int==2:
                            fb_dict[mod_int] = fb_dict[mod_int] & ( 0xFF - mask_int )
                if _feedback(fb_module.feedback,fb_mask.debugging):
                     sys.stderr.write(" feedback: mode %d on %d mask %d\n"%(
                         act_int,mod_int,mask_int))
    return r
