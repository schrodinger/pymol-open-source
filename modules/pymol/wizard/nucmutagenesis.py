import collections

from pymol.wizard import Wizard
import pymol.cmd
import pymol

SRC_SELE = "_mutate_sel"

FRAG_NAME = "_tmp_mut"

DEFAULT_MODE = "Adenine"
DEFAULT_REP = "lines"


class Status:
    NO_SELECTION = 0
    MUTAGENIZING = 1


class Nucmutagenesis(Wizard):
    _rep_name = {
        'lines': "Show Lines",
        'sticks': "Show Sticks",
        'spheres': "Show Spheres",
        'dots': "Show Dots",
    }

    _reps = list(_rep_name.keys())

    _auto_center_str = ["ON", "OFF"]

    # Nucleic Acid chi angle
    _all_dihedral_atoms = {
        "pur": ["O4'", "C1'", "N9", "C4"],
        "pyr": ["O4'", "C1'", "N1", "C2"]
    }

    _base_types = {
        "atp": "pur",
        "da": "pur",
        "a": "pur",
        "gtp": "pur",
        "dg": "pur",
        "g": "pur",
        "ctp": "pyr",
        "dc": "pyr",
        "c": "pyr",
        "ttp": "pyr",
        "dt": "pyr",
        "t": "pyr",  # non-canonical
        "utp": "pyr",
        "du": "pyr",  # non-canonical
        "u": "pyr"
    }

    _base3_to_1 = {"atp": "a", "ctp": "c", "gtp": "g", "ttp": "t", "utp": "u"}

    _base_atoms = [
        "N1", "HN1", "C2", "N2", "HN21", "HN22", "C2", "N3", "N4", "HN41",
        "HN42", "C4", "C5", "C5M", "HM51", "HM52", "HM53", "C6", "O6", "N6",
        "N7", "N8", "C8", "H8", "N9"
    ]

    _sugar_phos_atoms = [
        "C1'", "H1'", "O4'", "C2'", "H2'", "H2'1", "H2'2", "O2'", "HO2'",
        "C3'", "C4'", "H4'", "C5'", "H3'", "O3'", "HO3'", "H5'", "H5''",
        "H5'1", "O5", "PA", "O1A", "O2A", "HOA2", "O3A", "PB", "O1B", "O2B",
        "O3B", "HOB2", "O1G", "PG", "03G", "HOG3", "O2G", "HOG2", "O5'", "O3G",
        "P", "OP1", "OP2", "O1P", "O2P", "H1G", "H3G", "H1B", "H2A", "H5'2", "HB", "HA"
    ]

    _mode_labels = collections.OrderedDict([
        ('Adenine', 'ATP'),
        ('Cytosine', 'CTP'),
        ('Guanine', 'GTP'),
        ('Thymine', 'TTP'),
        ('Uracil', 'UTP'),
    ])

    def __init__(self, _self=pymol.cmd):  #PyMOL pattern
        Wizard.__init__(self, _self)  #PyMOL pattern
        if self.cmd.get_movie_length() > 0:
            raise pymol.wizarding.WizardError(
                'Mutagenesis Wizard cannot be used with Movie')

        self.cmd.unpick()

        self._stored = pymol.Scratch_Storage()
        self._space = {'stored': self._stored}

        self._status = Status.NO_SELECTION
        self._auto_center = "ON"
        self.mode = DEFAULT_MODE
        self.rep = DEFAULT_REP

        self.selection_mode = self.cmd.get_setting_int("mouse_selection_mode")
        self.cmd.set("mouse_selection_mode", 1)

        tmp_menu = [[2, 'Mutant', '']]
        for mode in self._mode_labels:
            tmp_menu.append(
                [1, mode, 'cmd.get_wizard().set_mode("' + mode + '")'])

        self.menu['mode'] = tmp_menu

        tmp_menu2 = [[2, 'Representation', '']]
        for rep in self._reps:
            tmp_menu2.append([
                1, self._rep_name[rep],
                'cmd.get_wizard().set_rep("' + rep + '")'
            ])
        self.menu['rep'] = tmp_menu2

        tmp_menu = [[2, 'Auto Center',
                     ''], [1, "ON", 'cmd.get_wizard().set_auto_center("ON")'],
                    [1, "OFF", 'cmd.get_wizard().set_auto_center("OFF")']]
        self.menu['auto_center'] = tmp_menu

    def get_prompt(self):
        '''
        Updates the upper-left prompt shown to guide user's next actions.
        Overrides method from Wizard class
        :return: Prompt text awaiting user's next wizard-related action
        :rtype: List (of size 1) of type string.
        '''

        self.prompt = None
        if self._status == Status.NO_SELECTION:
            self.prompt = ['Pick a nucleotide position...']
        elif self._status == Status.MUTAGENIZING:
            self.prompt = [
                'Select a nucleotide for %s or pick a new position...' %
                self.res_text
            ]
        return self.prompt

    def get_panel(self):
        '''
        Updates right-side panel to set/modify wizard-related actions
        Overrides method from Wizard class
        :return: List of widgets available for the wizard
        :rtype: List of [(widget-type), (widget text), (menu/command)] to be
                displayed on panel
        '''

        cmd = self.cmd
        if int(cmd.get_setting_int("mouse_selection_mode") != 1):
            cmd.set("mouse_selection_mode", 1)
        label = 'Mutate to ' + self.mode
        return [
            [1, 'Mutagenesis', ''],
            [3, label, 'mode'],
            [3, 'Auto Center: %s' % self._auto_center, 'auto_center'],
            [3, self._rep_name[self.rep], 'rep'],
            [2, 'Apply', 'cmd.get_wizard().apply()'],
            [2, 'Clear', 'cmd.get_wizard().clear()'],
            [2, 'Done', 'cmd.set_wizard()'],
        ]

    def do_pick(self, bondFlag):
        '''
        Pick callback for three-button editing during wizard
        Overrides method from Wizard class
        :param bondFlag: True if a bond is selected
        '''
        if bondFlag:
            self._error = "Please select an atom and not a bond."
            print(self._error)
        else:
            self.do_select('byres pk1')
        self.cmd.refresh_wizard()

    def do_select(self, selection):
        '''
        Selection callback during wizard
        Overrides method from Wizard class
        :param selection: A PyMOL selection
        :type selection: string
        '''
        cmd = self.cmd

        cmd.select(SRC_SELE, selection)
        cmd.unpick()
        try:
            self._do_mutation()
        except pymol.wizarding.WizardError as e:
            print(e)
        cmd.delete(selection)
        cmd.refresh_wizard()
        cmd.deselect()

    def _get_chi_dihedral_atoms(self, base_type_lower):
        '''
        Determines the atom types for the nucleic acid
        that composes the chi dihedral angle
        :param base_type_lower: base type in lower cases
        :type base_type_lower: string
        :return: a tuple of atom types which composes the
                chi dihedral
        :rtype: a tuple of four strings
        '''
        base_type = ""
        if base_type_lower in self._base_types:
            base_type = self._base_types[base_type_lower]
        else:
            if all(
                    self.cmd.count_atoms("(%s) & name %s" % (SRC_SELE,
                                                             name)) == 1
                    for name in self._all_dihedral_atoms["pur"]):
                base_type = "pur"
            else:
                base_type = "pyr"

        return tuple(self._all_dihedral_atoms[base_type])

    def _update_reps(self):
        '''
        Update atom representation of the candidate fragment
        '''
        self.cmd.show_as(self.rep + ' lines', '?' + FRAG_NAME)

    def _transfer_dihedral(self):
        '''
        Calculates the dihedral angle from the original residue
        and applies it to the candidate fragment
        '''
        chi_dihedral_src = self.cmd.get_dihedral(
            "(%s & name %s)" % (SRC_SELE, self._src_O_atom_name),
            "(%s & name %s)" % (SRC_SELE, self._src_Cp_atom_name),
            "(%s & name %s)" % (SRC_SELE, self._src_N_atom_name),
            "(%s & name %s)" % (SRC_SELE, self._src_C_atom_name))

        self.cmd.set_dihedral(
            "(%s & name %s)" % (FRAG_NAME, self._frag_O_atom_name),
            "(%s & name %s)" % (FRAG_NAME, self._frag_Cp_atom_name),
            "(%s & name %s)" % (FRAG_NAME, self._frag_N_atom_name),
            "(%s & name %s)" % (FRAG_NAME, self._frag_C_atom_name),
            chi_dihedral_src)

    def _do_mutation(self):
        '''
        After selection, propose what the fragment will look like where
        applied.
        Exception: _do_mutation should be wrapped in a try-exception block
                   when called
        '''
        cmd = self.cmd

        cmd.delete(FRAG_NAME)
        cmd.refresh_wizard()

        if any(
                cmd.count_atoms("(%s) & name %s" % (SRC_SELE, name)) != 1
                for name in ("C1'", "C2'", "C3'", "C4'", "O4'")):
            self.clear()
            raise pymol.wizarding.WizardError(
                'Improper selection of nucleic acid.')

        frag_type_lower = self._mode_labels[self.mode].lower()
        cmd.fragment(frag_type_lower, FRAG_NAME, origin=0)
        self._update_reps()
        cmd.color("white", "%s & elem C" % FRAG_NAME)
        cmd.iterate(
            "first (%s)" % SRC_SELE,
            'stored.name = model+"/"+segi+"/"+chain+"/"+resn+"`"+resi',
            space=self._space)
        self.res_text = self._space['stored'].name

        cmd.alter(
            "?%s & name C1'" % SRC_SELE,
            "stored.identifiers = (segi, chain, resi, ss, color)",
            space=self._space)
        self._status = Status.MUTAGENIZING
        self.get_prompt()

        cmd.iterate(
            "(%s & name C1')" % SRC_SELE,
            "stored.resn = resn",
            space=self._space)
        src_type_lower = self._space['stored'].resn.lower()
        self._src_O_atom_name, self._src_Cp_atom_name, self._src_N_atom_name, \
        self._src_C_atom_name = self._get_chi_dihedral_atoms(src_type_lower)

        self._frag_O_atom_name, self._frag_Cp_atom_name, self._frag_N_atom_name, \
        self._frag_C_atom_name = self._get_chi_dihedral_atoms(frag_type_lower)

        cmd.pair_fit("(%s & name %s)" % (FRAG_NAME, self._frag_O_atom_name),
                     "(%s & name %s)" % (SRC_SELE, self._src_O_atom_name),
                     "(%s & name %s)" % (FRAG_NAME, self._frag_Cp_atom_name),
                     "(%s & name %s)" % (SRC_SELE, self._src_Cp_atom_name),
                     "(%s & name %s)" % (FRAG_NAME, self._frag_N_atom_name),
                     "(%s & name %s)" % (SRC_SELE, self._src_N_atom_name))

        self._transfer_dihedral()

        for a in self._sugar_phos_atoms:
            cmd.remove("(%s & name %s)" % (FRAG_NAME, a))

        if self._auto_center == "ON":
            cmd.center(FRAG_NAME, animate=1)

    def clear(self):
        '''
        Clear selections and objects created by this wizard
        in PyMOL.
        Overrides method from Wizard class
        '''
        cmd = self.cmd
        self._status = Status.NO_SELECTION
        self.mode = DEFAULT_MODE
        self.rep = DEFAULT_REP
        cmd.delete(FRAG_NAME)
        cmd.refresh_wizard()

    def set_mode(self, mode):
        '''
        Modifies the wizard's mode aka chooses which fragment becomes
        the candidate nucleic acid fragment
        :param mode: Wizard-specific mode. Here, modes are any of the
                    nucleic acid fragments
        :type mode: string
        '''
        cmd = self.cmd

        #avoid case-sensitivity
        mode_fixed = mode.title()
        if mode_fixed in self._mode_labels:
            self.mode = mode_fixed
        else:
            print('Improper Nucleic Acid')
        if self._status == Status.MUTAGENIZING:
            try:
                self._do_mutation()
            except pymol.wizarding.WizardError as e:
                print(e)
        cmd.refresh_wizard()

    def set_rep(self, rep):
        '''
        Sets the atom representation of the candidate fragment
        Overrides method from Wizard class
        :param rep: atom representation
        :rtype rep: string
        '''
        cmd = self.cmd
        if self._status != Status.MUTAGENIZING:
            return
        if rep in self._reps:
            self.rep = rep
        self._update_reps()
        cmd.refresh_wizard()

    def set_auto_center(self, auto_center):
        '''
        Sets the flag for automatic centering upon prospective mutation
        :param auto_center: string flag ("ON" or "OFF")
        :rtype auto_center: string
        '''
        cmd = self.cmd
        if auto_center not in self._auto_center_str:
            print("Improper Auto Center setting. 'ON' or 'OFF' accepted only")
        self._auto_center = auto_center
        cmd.refresh_wizard()

    def cleanup(self):
        '''
        Reverts the wizard to its original state
        Overrides method from Wizard class
        '''
        cmd = self.cmd
        cmd.set("mouse_selection_mode",
                self.selection_mode)  # restore selection mode
        self.clear()

    def apply(self):
        '''
        Permanently applies the mutation upon the selected residue
        Overrides method from Wizard class
        '''
        cmd = self.cmd
        if self._status == Status.NO_SELECTION:
            return

        #Remove all atoms that are not in the sugar/phosphate group
        #Needed for non-canonical bases
        cmd.select("_tmp_sele_invert", "(none)")
        for a in self._sugar_phos_atoms:
            cmd.select("_tmp_sele_invert",
                       "_tmp_sele_invert | (%s & name %s)" % (SRC_SELE, a))

        cmd.select("_tmp_sele_invert",
                   "%s and not _tmp_sele_invert" % SRC_SELE)

        cmd.remove("_tmp_sele_invert")

        try:
            new_name = cmd.get_object_list(SRC_SELE)[0]
        except IndexError:
            print(" Mutagenesis: object not found.")
            return

        frag_name_three = self._mode_labels[self.mode]
        frag_name_one = self._base3_to_1[frag_name_three.lower()]

        # FUSE
        cmd.fuse(
            "%s & name %s" % (FRAG_NAME, self._frag_N_atom_name),
            "/%s/%s/%s/%s & name %s" %
            (new_name, self._stored.identifiers[0],
             self._stored.identifiers[1], self._stored.identifiers[2],
             self._src_Cp_atom_name), 1)

        #Check to see if DNA
        dnaPrefix = ''
        if cmd.count_atoms("%s & name O2'" % SRC_SELE) == 0:
            dnaPrefix = 'D'

        cmd.alter(
            "/%s/%s/%s/%s" %
            (new_name, self._stored.identifiers[0],
             self._stored.identifiers[1], self._stored.identifiers[2]),
            "(resn) = '%s%s'" % (dnaPrefix, frag_name_one.upper()),
            space=self._space)

        cmd.unpick()
        self.clear()
