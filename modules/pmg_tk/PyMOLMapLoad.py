# TODO
# * default settings for color and rep
# * make the final viewing step a function
# *   setup_map(name,levels=,colors=,reps=)

import pymol
from pymol import headering
import Pmw

class PyMOLMapLoad:

    def __init__(self,parent,app,f):
        self._parent   = parent
        self._app      = app
        self._fileName = f
        self._fileData = None

        # model state
        self._amplitudes = None
        self._phases     = None
        self._weights    = None
        self._min_res    = None
        self._max_res    = None
        self._fofc       = None
        self._name_prefix= None

        # reflection file header data
        if f[-3:] in ("MTZ", "mtz"):
            self._fileData = headering.MTZHeader(f)
        elif f[-3:] in ("CIF", "cif"):
            self._fileData = headering.CIFHeader(f)
        elif f[-3:] in ("CNS", "cns", "hkl", "HKL"):
            self._fileData = headering.CNSHeader(f)

    def pack_and_show(self):
        self.pack()
        return self.show()

    def pack(self):
        # MAIN DIALOG
        self._d = Pmw.Dialog(self._parent,
                             buttons = ("OK", "Cancel", "Help"),
                             defaultbutton = "OK",
                             title = "PyMOL Map Generation",
                             command = self.run)
        self._d.geometry("+%d+%d" % (self._app.winfo_reqwidth(), self._app.winfo_reqheight()))
        self._d.withdraw()
        self._d.protocol('WM_DELETE_WINDOW', self.quit)

        #
        # COLUMN LABEL GROUP
        #
        self._col_gp = Pmw.Group(self._d.interior(),
                                 tag_text="Column Labels",)
        self._col_gp.pack(fill='x', expand='yes')

        defaultListHeight = 125
        FCols = []
        FCols.extend(self._fileData.getColumnsOfType("F"))
        FCols.extend(self._fileData.getColumnsOfType("G"))
        if not len(FCols): FCols = [ "" ]
        self._ampl_chooser = Pmw.ComboBox(self._col_gp.interior(),
                                          label_text = "Amplitudes",
                                          labelpos = "nw",
                                          selectioncommand = self.set_amplitudes, 
                                          scrolledlist_items = FCols,             
                                          dropdown = 1,
                                          listheight=defaultListHeight,
                                          sticky='ew')
        self._ampl_chooser.pack(fill='both',expand=1,padx=7,pady=4)
        _FC, _PC, _looksLike = self._fileData.guessCols("FoFc")
        _2FC, _2PC, _looksLike = self._fileData.guessCols("2FoFc")
        # be nice and choose the most appropriate col
        if _2FC!=None:
            if _2FC in FCols:
                self._ampl_chooser.selectitem(_2FC)
        elif _FC!=None:
            if _FC in FCols:
                self._ampl_chooser.selectitem(_FC)
        else:
            self._ampl_chooser.selectitem(FCols[0])

        PCols = []
        PCols.extend(self._fileData.getColumnsOfType("P"))
        if not len(PCols): PCols = [ "" ]
        self._phase_chooser = Pmw.ComboBox(self._col_gp.interior(),
                                          label_text = "Phases",
                                          labelpos = "nw",
                                          selectioncommand = self.set_phases, 
                                          scrolledlist_items = PCols,             
                                          dropdown = 1,
                                          listheight=defaultListHeight)
        self._phase_chooser.pack(fill='both', expand=1,padx=7,pady=4)
        # be nice and choose the most appropriate col
        if _2PC!=None:
            if _2PC in PCols:
                self._phase_chooser.selectitem(PCols.index(_2PC))
        elif _PC!=None:
            if _PC in PCols:
                self._phase_chooser.selectitem(PCols.index(_PC))
        else:
            self._phase_chooser.selectitem(PCols[0])

        WCols = [ "None", ]
        WCols.extend(self._fileData.getColumnsOfType("W"))
        WCols.extend(self._fileData.getColumnsOfType("Q"))
        self._wt_chooser = Pmw.ComboBox(self._col_gp.interior(),
                                          label_text = "Weights",
                                          labelpos = "nw",
                                          selectioncommand = self.set_weights,
                                          scrolledlist_items = WCols,             
                                          dropdown = 1,
                                          listheight=defaultListHeight)
        self._wt_chooser.pack(fill='both', expand=1,padx=7,pady=4)
        self._wt_chooser.selectitem("None")

        #
        # INPUT OPTIONS GROUP
        #
        self._input_gp = Pmw.Group(self._d.interior(),
                                   tag_text="Input Options",)
        self._input_gp.pack(fill='both', expand='yes')

        if self._fileData.reso_min!=None:
            default_min_res = float("%3.5f"%float(self._fileData.reso_min))
        else:
            default_min_res = ""
        if self._fileData.reso_max!=None:
            default_max_res = float("%3.5f"%float(self._fileData.reso_max))
        else:
            default_max_res = ""
        self._min_res_fld = Pmw.EntryField(self._input_gp.interior(),
                                           labelpos="wn",
                                           label_text="Min. Resolution",
                                           value = default_min_res,
                                           validate = { "validator" : 'real' },
                                           entry_width=7,
                                           modifiedcommand=self.set_min_res,
                                           command = self.set_min_res)
        self._min_res_fld.grid(row=1,column=0,rowspan=2,sticky='ew',pady=4)

        self._max_res_fld = Pmw.EntryField(self._input_gp.interior(),
                                           labelpos="wn",
                                           label_text = "Max Resolution",
                                           value = default_max_res,
                                           validate = { "validator" : 'real' },
                                           entry_width=7,
                                           modifiedcommand=self.set_max_res,
                                           command = self.set_max_res)
        self._max_res_fld.grid(row=1,column=1,rowspan=2,sticky='ew',pady=4)


        #
        # MAP OPTIONS GROUP
        #
        self._options_gp = Pmw.Group(self._d.interior(),
                                     tag_text="Map Options",)
        self._options_gp.pack(fill='x', expand='yes')                                

        self._name_prefix_fld = Pmw.EntryField(self._options_gp.interior(),
                                               labelpos="wn",
                                               label_text = "New Map Name Prefix",
                                               value = "",
                                               validate = { "validator" : 'alphanumeric' },
                                               entry_width=20,
                                               modifiedcommand=self.set_name_prefix,
                                               command = self.set_name_prefix)
        self._name_prefix_fld.pack(fill="x", expand=0, anchor='w')

        self._fofc_chooser = Pmw.RadioSelect(self._options_gp.interior(),
                                             command = self.set_fofc,
                                             buttontype="checkbutton",)
        self._fofc_chooser.add("FoFc")
        self._fofc_chooser.pack(fill="none", expand=0, anchor="w")

    def show(self):
        self._d.show()
    def quit(self):
        if __name__=="__main__":
            # TODO--remove me; use for development only!
            self._parent.destroy()
        else:
            # TODO -- use only this in release
            self._d.destroy()
        


    # UI SETTERS
    def set_amplitudes(self,arg):
        self._amplitudes = arg
    def set_phases(self,arg):
        self._phases = arg
    def set_weights(self,arg):
        self._weights = arg
    def set_min_res(self):
        self._min_res = self._min_res_fld.getvalue()
    def set_max_res(self):
        self._max_res = self._max_res_fld.getvalue()
    def set_fofc(self,arg,state):
        self._fofc = state                                      
    def set_name_prefix(self):
       self._name_prefix = self._name_prefix_fld.getvalue()

    def update_state(self):
        # grab all values
        self._amplitudes = self._ampl_chooser.get()
        self._phases     = self._phase_chooser.get()
        self._weights    = self._wt_chooser.get()
        self._min_res    = self._min_res_fld.getvalue()
        self._max_res    = self._max_res_fld.getvalue()
        self._fofc       = len(self._fofc_chooser.getvalue())>0
        self._name_prefix= self._name_prefix_fld.getvalue()
        
    def report_state(self):
        print("Here is the state of the box")
        print("Amplitudes:\t%s" % self._amplitudes)
        print("Phases    :\t%s" % self._phases)
        print("Weights   :\t%s" % self._weights)
        print("Min Res   :\t%s" % self._min_res)
        print("Max Res   :\t%s" % self._max_res)
        print("FoFc      :\t%s" % str(self._fofc))
        print("Name Prefix :\t'%s'" % self._name_prefix)

    def show_help(self,msg=None,title=None):
        # TODO -- CHANGE THE HELP TEXT
        if msg==None:
            helpText = pymol.cmd.map_generate.__doc__
        else:
            helpText = msg

        if title==None:
            title="PyMOL Map Loading Help"

        h = Pmw.TextDialog(self._parent,
                           title=title,)
        h.insert("end", helpText)
        h.configure(text_state='disabled')


    def run(self,action):
        if action=="OK":
            self.update_state()
            #self.report_state()

            if self._name_prefix==None or self._name_prefix=="":
                # grep the dataset name from amplitudes
                if '/' in self._amplitudes:
                    pfx = self._amplitudes.split('/')
                    if len(pfx)>=2:
                        pfx = pfx[1]
                else:
                    pfx = self._amplitudes
            else:
                pfx = self._name_prefix

            # to ensure a clean name
            pfx = pymol.cmd.get_unused_name(pfx)

            if not len(self._amplitudes):
                missing_ampl = """
To synthesize a map from reflection data you need to specify at
leastone column for amplitudes and one column for phases. The
amplitudes column name was blank, and therefore PyMOL cannot create
the map.  Please select an amplitude column name from the file and try
again.
               """
                self.show_help(missing_ampl,"Missing Amplitudes Column Name")
                return None
            
            if not len(self._phases):
                missing_phases = """
To synthesize a map from reflection data you need to specify at least
one column for amplitudes and one column for phases. The phases column
name was blank, and therefore PyMOL cannot create the map.  Please
select an amplitude column name from the file and try again.
               """
                self.show_help(missing_phases, "Missing Phases Column Name")
                return None
            
            try:
                r = pymol.cmd.map_generate(pfx, self._fileName,
                                       self._amplitudes, self._phases, self._weights,
                                       self._min_res, self._max_res, 1, 1)
            except pymol.CmdException as e:
                print(e)
                return None

            if r==None or r=="None" or r=="":
                print(" MapLoad-Error: PyMOL could not load the MTZ file '%s' due to an unspecified error." % self._fileName)
                print(" MapLoad-Error: This typically occurs with bad data or blank column names. Please try again")
                print(" MapLoad-Error: or contact 'help@schrodinger.com' for more information.")
                return None
            skin     = pymol._ext_gui.skin
            try:
                pymol.cmd.set("suspend_updates", 1)
                if self._fofc:
                    toShow = pymol.cmd.get_setting_text("default_fofc_map_rep")
                    if toShow=="isosurface":
                        pymol.cmd.isosurface(pymol.cmd.get_unused_name(r+"-srf"),
                                             pfx, level=1.0)
                    elif toShow=="isomesh":
                        meshName=pymol.cmd.get_unused_name(r+"-msh3")
                        pymol.cmd.isomesh(meshName, pfx, level=3.0)
                        pymol.cmd.color("green", meshName)
            
                        meshName=pymol.cmd.get_unused_name(r+"-msh-3")
                        pymol.cmd.isomesh(meshName, pfx, level=-3.0)
                        pymol.cmd.color("red", meshName)
                    else:
                        # setup volume view
                        volName = pymol.cmd.get_unused_name(r+"-vol")
                        pymol.cmd.volume(volName, pfx, "fofc")
                        # if you don't do this, PyMOL will crash
                        # when it tries to load the panel
                else:
                    toShow = pymol.cmd.get_setting_text("default_2fofc_map_rep")
                    if toShow=="isosurface":
                        surfName=pymol.cmd.get_unused_name(r+"-srf")
                        pymol.cmd.isosurface(surfName, pfx, level=1.0)
                        pymol.cmd.color("blue", surfName)
                    elif toShow=="isomesh":
                        meshName=pymol.cmd.get_unused_name(r+"-msh")
                        pymol.cmd.isomesh(meshName, pfx, level=1.0)
                        pymol.cmd.color("blue", meshName)
                    else:
                        # setup volume view
                        volName = pymol.cmd.get_unused_name(r+"-vol")
                        pymol.cmd.volume(volName, pfx, "2fofc")
                        # if you don't do this, PyMOL will crash
                        # when it tries to load the panel
                        
            except:
                pass
            finally:
                pymol.cmd.set("suspend_updates", 0)

            if r!=None:
                # setting?
                if pymol.cmd.get_setting_boolean("autoclose_dialogs"):
                    self.quit()

        elif action=="Cancel":
            self.quit()
        elif action=="Help":
            self.show_help()
