from pymol.wizard import Wizard
from pymol import cmd
import pymol
import traceback

sele_prefix = "_pf_s_"
sele_prefix_len = len(sele_prefix)

dist_prefix = "_pf_d_"

indi_sele = "_indicate_pf"

class Pair_fit(Wizard):

    def __init__(self,_self=cmd):

        Wizard.__init__(self,_self)

        self.memory = 0
        self.n_pair = 0
        self.status = 0 # 0 no atoms selections, 1 atom selected
        self.message = None

        self.selection_mode = cmd.get_setting_int("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect() # disable the active selection (if any)

    def get_panel(self):
        return [
            [ 1, 'Pair Fitting',''],
            [ 2, 'Fit %d Pairs'%self.n_pair,'cmd.get_wizard().fit()'],
            [ 2, 'Delete Last Pair','cmd.get_wizard().remove_last()'],
            [ 2, 'Redraw','cmd.get_wizard().update_dashes()'],
            [ 2, 'Clear','cmd.get_wizard().clear()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def cleanup(self):
        self.clear()
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode

    def clear(self):
        cmd.delete(sele_prefix+"*")
        cmd.delete(dist_prefix+"*")
        cmd.delete(indi_sele)
        lst = cmd.get_names('selections')
        self.n_pair = 0
        self.status = 0
        self.message = None
        cmd.unpick()
        cmd.refresh_wizard()

    def get_prompt(self):
        self.prompt = None
        if self.status==0:
            self.prompt = [ 'Pick the mobile atom...']
        elif self.status==1:
            self.prompt = [ 'Pick the target atom...' ]
        if self.message is not None:
            self.prompt.append(self.message)
        return self.prompt

    def set_status(self,status):
        self.status = status
        cmd.refresh_wizard()

    def get_sele_list(self,mode='all'):
        lst = cmd.get_names('selections')
        lst = [x for x in lst if x[0:sele_prefix_len]==sele_prefix]
        lst.sort()
        if mode == 'mobile': # mobile
            lst=[x for x in lst if x[-1:]=='b']
        elif mode == 'target': # target
            lst=[x for x in lst if x[-1:]=='a']
        return lst

    def fit(self):
        # build up the pair-wise list of selections
        cmd.delete(dist_prefix+"*")
        lst = self.get_sele_list()
        c = 0
        args = []
        while 1:
            if not len(lst): break
            a = lst.pop()
            if not len(lst): break
            b = lst.pop()
            args.append(a)
            args.append(b)
        # do the fit
        if len(args):
            cmd.push_undo(args[0])
            dist = cmd.pair_fit(*args)
            self.message = "RMS over %d pairs = %5.3f"%(self.n_pair,dist)
            cmd.refresh_wizard()
        self.update_dashes()

    def remove_last(self):
        # build up the pair-wise list of selections
        cmd.delete(dist_prefix+"*")
        lst = self.get_sele_list()
        if len(lst):
            cmd.delete(lst.pop())
            if len(lst):
                cmd.delete(lst.pop())
            self.n_pair = self.n_pair - 1
        self.update_dashes()
        self.status=0
        cmd.refresh_wizard()

    def update_dashes(self):
        cmd.delete(dist_prefix+"*")
        lst = self.get_sele_list()
        c = 0
        while 1:
            if not len(lst): break
            a = lst.pop()
            if not len(lst): break
            b = lst.pop()
            name = dist_prefix+str(c)
            cmd.dist(name,a,b,width=7,length=0.05,gap=0.05)
            cmd.hide('label',name)
            cmd.enable(name)
            c = c + 1

    def check_same_object(self,lst,sele):
        if not len(lst):
            return 1
        else:
            if cmd.count_atoms("((byobj %s) and %s)"%(lst[0],sele),quiet=1):
                return 1
        return 0

    def check_different_object(self,lst,sele):
        if not len(lst):
            return 1
        else:
            if not cmd.count_atoms("((byobj %s) and %s)"%(lst[0],sele),quiet=1):
                return 1
        return 0

    def do_select(self,name): # map selects into picks
        cmd.unpick()
        try:
            cmd.edit(name + " and not " + sele_prefix + "*") # note, using new object name wildcards
            cmd.delete(name)
            self.do_pick(0)
        except pymol.CmdException:
            traceback.print_exc()
            pass

    def do_pick(self,bondFlag):
        if bondFlag:
            self.message = "Error: please select an atom, not a bond."
            print(self.message)
        else:
            if self.status==0:
                lst = self.get_sele_list(mode='mobile')
                if not self.check_same_object(lst,"(pk1)"):
                    self.message = "Error: must select an atom in the same object as before."
                    print(self.message)
                else:
                    name = sele_prefix + "%02db"%self.n_pair # mobile end in 'b'
                    cmd.select(name,"(pk1)")
                    cmd.unpick()
                    cmd.select(indi_sele,name)
                    cmd.enable(indi_sele)
                    self.status = 1
                    self.message = None
            elif self.status==1:
                lst = self.get_sele_list(mode='target')
                if not self.check_same_object(lst,"(pk1)"):
                    self.message = "Error: must select an atom in the same object as before."
                    print(self.message)
                else:
                    lst = self.get_sele_list(mode='mobile')
                    if not self.check_different_object(lst,"(pk1)"):
                        self.message = "Error: target atom must be in a distinct object."
                        print(self.message)
                    else:
                        name = sele_prefix + "%02da"%self.n_pair # target end in 'a'
                        cmd.select(name,"(pk1)")
                        cmd.unpick()
                        cmd.select(indi_sele,name)
                        cmd.enable(indi_sele)
                        self.n_pair = self.n_pair + 1
                        self.status = 0
                        self.update_dashes()

        cmd.refresh_wizard()
