# filter wizard
# no-frills tool for quickly filtering docked compounds, etc.

import os,sys
from pymol.wizard import Wizard
from pymol import cmd
import traceback
# global dictionary for saving result on a per-object basis

static_dict = {}

# last/current object being filtered

default_object = None

# browing mode

default_browse = 1

accept_str = "Accept"
defer_str = "Defer"
reject_str = "Reject"

# class definition (class name must match wizard name with cap)

class Filter(Wizard):

    def migrate_session(self, version):
        if version >= self.cmd.get_version()[2]:
            return

        # remap old "title" state identifiers to new identifiers
        # (changed in 1.7.0.0)

        for object, sdo in self.dict.items():
            if not sdo:
                continue

            # build title -> ident mapping
            tota = self.cmd.count_states(object)
            title2ident = dict(
                    (self.cmd.get_title(object, state),
                        self.get_ident(object, state))
                    for state in range(1, tota + 1))

            # build remapped sdo
            try:
                new_sdo = {}
                for title, value in sdo.items():
                    new_sdo[title2ident[title]] = value
                sdo.clear()
                sdo.update(new_sdo)
            except KeyError:
                # lookup by title failed, we can assume that this instance
                # already uses new identifiers
                return

        self.load_state_dict()

    def __init__(self,_self=cmd):

        # initialize parent class

        Wizard.__init__(self,_self)

        self.update_object_menu()

        # restore previous state from global storage

        self.dict = static_dict
        self.object = default_object if default_object in self.avail_objects else None
        self.browse = default_browse
        self.state_dict = {}

        # if we don't have a current object, choose the first multi-state object

        if not self.object and self.avail_objects:
            self.object = self.avail_objects[0]

        # menu for

        self.menu['browse'] = [
            [2, 'Browse Mode',''],
            [1, 'Browse All','cmd.get_wizard().set_browse(1)'],
            [1, 'Browse Accepted','cmd.get_wizard().set_browse(2)'],
            [1, 'Browse Rejected','cmd.get_wizard().set_browse(3)'],
            [1, 'Browse Deferred','cmd.get_wizard().set_browse(4)'],
            [1, 'Browse Remaining','cmd.get_wizard().set_browse(5)'],
            ]

        self.menu['create'] = [
            [2, 'Create Filtered Object', ''],
            [1, 'Accepted','cmd.get_wizard().create_object("Accept")'],
            [1, 'Rejected','cmd.get_wizard().create_object("Reject")'],
            [1, 'Deferred','cmd.get_wizard().create_object("Defer")'],
            ]

        self.count_object()
        self.load_state_dict()
        self.update_object_menu()
        cmd.set_key('F1',lambda s=self:s.accept())
        cmd.set_key('F2',lambda s=self:s.reject())
        cmd.set_key('F3',lambda s=self:s.defer())
        cmd.set_key('right',lambda s=self:s.forward())
        cmd.set_key('left',lambda s=self:s.backward())

    def do_select(self,name):
        try:
            obj_name = cmd.index('first ?' + name)[0][0]
            self.set_object(obj_name)
        except:
            pass
        self.cmd.deselect()

    def do_pick(self,bondFlag):
        self.do_select('pk1')
        self.cmd.unpick()

    def do_state(self, state):
        cmd.refresh_wizard()

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select + Wizard.event_mask_state

    def update_object_menu(self):

        # find objects with > 1 state

        self.avail_objects = []

        for a in cmd.get_names('objects'):
            if cmd.get_type(a)=='object:molecule':
                if cmd.count_states(a)>1:
                    self.avail_objects.append(a)

        # now create the object menu

        self.menu['object'] = [[2,'Select Object','']]
        for a in self.avail_objects:
            self.menu['object'].append([ 1,a,'cmd.get_wizard().set_object("%s")'%(a) ])
        self.menu['object'].append([ 1,'None','cmd.get_wizard().set_object(None)'])

    def set_browse(self,browse):
        # allow user to focus on only a subset of the compounds
        self.browse = browse
        if self.browse == 1:
            print(" Filter: Browsing all compounds.")
            cmd.mset() # all states visible
        elif self.object is None:
            print(" Filter-Error: please choose an object first")
        else:
            self.check_object_dict()
            if self.browse == 2:
                print(" Filter: Browsing accepted compounds.")
                target = accept_str
            elif self.browse == 3:
                print(" Filter: Browsing rejected compounds.")
                target = reject_str
            elif self.browse == 4:
                print(" Filter: Browsing deferred compounds.")
                target = defer_str
            lst = []
            sd = self.state_dict
            sdo = self.dict[self.object]
            if self.browse<5:
                for a in list(sdo.keys()):
                    if sdo[a]==target:
                        lst.append(sd[a])
            else:
                print(" Filter: Browsing remaining compounds")
                for a in sd.keys():
                    if a not in sdo:
                        lst.append(sd[a])
            lst.sort()
            if len(lst)==0:
                print(" Filter-Error: No matching compounds.")
            cmd.mset(' '.join(map(str,lst)))
            cmd.rewind()
        cmd.refresh_wizard()

    def check_object_dict(self):
        # make sure we have a valid entry for this object in our dictionary

        if self.object not in self.dict:
            self.dict[self.object]={} # create dictionary to store results

    def adjust(self,decision,inc):
        # utility routine to increment/decrement counters
        if decision == accept_str:
            self.acce = self.acce + inc
        elif decision == reject_str:
            self.reje = self.reje + inc
        elif decision == defer_str:
            self.defe = self.defe + inc

    def load_state_dict(self):
        # establish relationship between names and states
        # ASSUMPTION: identifiers will be unique
        self.state_dict = {}
        sd = self.state_dict
        so = self.object
        if so is not None:
            cnt = cmd.count_states(so)
            for a in range(1,cnt+1):
                sd[self.get_ident(so,a)] = a

    def count_object(self):
        # record how many molecular are in an object, etc.
        self.check_object_dict()
        if self.object is not None:
            self.acce = 0
            self.reje = 0
            self.defe = 0
            self.togo = 0
            self.tota = cmd.count_states(self.object)
            sdo=self.dict[self.object]
            self.togo = self.tota-len(sdo)
            for a in list(sdo.keys()):
                dec = sdo[a]
                self.adjust(dec,1)

    def set_object(self,obj_name):
        self.object = obj_name
        self.count_object()
        self.load_state_dict()
        cmd.refresh_wizard()

    def get_panel(self):

        # returns Wizard panel for PyMOL to display

        # 1 = title/text
        # 2 = button
        # 3 = pop-up menu

        self.update_object_menu()
        if self.object is not None:
            save_str = 'Save %s.txt'%self.object
        else:
            save_str = ""
        return [
            [ 1, 'Filtering Wizard',''],
            [ 3, self.menu['browse'][self.browse][1], 'browse' ],
            [ 3, 'Object: %s' % (self.object), 'object' ],
            [ 2, 'Accept (F1)','cmd.get_wizard().accept()'],
            [ 2, 'Reject (F2)','cmd.get_wizard().reject()'],
            [ 2, 'Defer (F3)','cmd.get_wizard().defer()'],
            [ 2, 'Forward (->)','cmd.get_wizard().forward()'],
            [ 2, 'Back (<-)','cmd.get_wizard().backward()'],
            [ 3, 'Create Filtered Object', 'create'],
            [ 2, save_str,'cmd.get_wizard().save()'],
            [ 2, 'Refresh','cmd.refresh_wizard()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def get_ident(self, object, state):
        return '%d/%d %s' % (state, self.tota,
                self.cmd.get_title(self.object, state))

    def get_prompt(self):

        # returns text prompt

        self.prompt = None
        if self.object is None:
            self.prompt = [ 'Please select a multi-state object...' ]
        else:
            self.prompt = [ '%s: %d accepted, %d rejected, %d deferred, %d remaining'%(
                self.object,self.acce,self.reje,self.defe,self.togo) ]
            state = cmd.get_object_state(self.object)
            ident = self.get_ident(self.object,state)
            sdo=self.dict[self.object]
            if ident in sdo:
                self.prompt.append('%s: %s'%(ident,sdo[ident]))
            else:
                self.prompt.append('%s?'%(ident))
        return self.prompt

    def count(self,entry,str):
        # keep track of how many compounds are in which category

        self.check_object_dict()
        sdo = self.dict[self.object]
        if entry in sdo:
            self.adjust(sdo[entry],-1)
        else:
            self.togo = self.togo - 1
        sdo[entry] = str
        self.adjust(sdo[entry],1)

    def accept(self):
        # accept compound and advance
        if self.object is None:
            print(" Filter-Error: Please choose an object first")
        else:
            state = cmd.get_object_state(self.object)
            ident = self.get_ident(self.object,state)
            print(" Filter: Accepting '%s'"%ident)
            self.count(ident,accept_str)
        cmd.forward()
        cmd.refresh_wizard()

    def reject(self):
        # reject compound and advance
        if self.object is None:
            print(" Filter-Error: Please choose an object first")
        else:
            state = cmd.get_object_state(self.object)
            ident = self.get_ident(self.object,state)
            print(" Filter: Rejecting '%s'"%ident)
            self.check_object_dict()
            self.count(ident,reject_str)
        cmd.forward()
        cmd.refresh_wizard()

    def defer(self):
        # defer compound and advance
        if self.object is None:
            print(" Filter-Error: Please choose an object first")
        else:
            state = cmd.get_object_state(self.object)
            ident = self.get_ident(self.object,state)
            print(" Filter: Deferring '%s'"%ident)
            self.check_object_dict()
            self.count(ident,defer_str)
        cmd.forward()
        cmd.refresh_wizard()

    def forward(self):
        # go forward and update information
        cmd.forward()
        cmd.refresh_wizard()

    def backward(self):
        # go backward and update information
        cmd.backward()
        cmd.refresh_wizard()

    def create_object(self, what='Accept'):
        if not self.object:
            print(" Filter-Error: Please choose an object first")
            return
        name = self.cmd.get_unused_name(self.object + '_' + what, 0)
        sdo = self.dict[self.object]
        lst = [self.state_dict[ident] for (ident, w) in sdo.items() if w == what]
        for state in sorted(lst):
            self.cmd.create(name, self.object, state, -1)

    def save(self):
        # write compounds to a file
        if self.object is None:
            print(" Filter-Error: please choose an object first")
        else:
            self.check_object_dict()
            fname = self.object+".txt"
            try:
                f=open(fname,'w')
                f.close()
            except:
                print(" Filter-Warning: '"+fname+"' in current directory is not writable.")
                print(" Filter-Warning: attempting to write in home directory.")
                fname = cmd.exp_path(os.path.join('~', fname))
            try:
                f=open(fname,'w')
                sd = self.state_dict
                sdo = self.dict[self.object]
                f.write('Object\t"%s"\n'%(self.object))
                f.write('Total\t%d\nAccepted\t%d\nRejected\t%d\nDeferred\t%d\nRemaining\t%d\n\n'%(
                    self.tota,
                    self.acce,
                    self.reje,
                    self.defe,
                    self.togo))
                # sort output in order of states
                lst = []
                for a in sd.keys():
                    lst.append((sd[a],a))
                lst.sort()
                # write list with decisions
                for a in lst:
                    if a[1] in sdo:
                        f.write('%d\t"%s"\t"%s"\n'%(a[0],a[1],sdo[a[1]]))
                    else:
                        f.write('%d\t"%s"\t"?"\n'%(a[0],a[1]))
                f.close()
                print(" Filter: Wrote '%s'."%fname)
            except:
                traceback.print_exc()
                print(" Filter-Error: Unable to write '%s'."%fname)

    def cleanup(self):
        # save current state in global vars...
        global default_object,default_browse,static_dict
        default_object = self.object
        default_browse = self.browse
        static_dict = self.dict
        # restore key actions
        cmd.set_key('F1',None)
        cmd.set_key('F2',None)
        cmd.set_key('F3',None)
        cmd.set_key('right',cmd.forward)
        cmd.set_key('left',cmd.backward)
