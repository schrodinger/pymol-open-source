'''
unit tests for pymol.exporting
'''

import os
import sys
import tempfile
import Image
import unittest

import pymol.exporting
from pymol import cmd, testing, stored

class TestExporting(testing.PyMOLTestCase):

    def testCache(self):
        for action in pymol.exporting.cache_action_dict:
            cmd.cache(action)

    def testCopyImage(self):
        cmd.copy_image
        self.skipTest("TODO")

    def testGetFastastr(self):
        seq, name = 'ACD', 'm1'
        cmd.fab(seq, name)
        s = cmd.get_fastastr()
        lines = s.split()
        self.assertTrue(lines[0] in (
            '>' + name,
            '>' + name + '_',
        ))
        self.assertEqual(lines[1:], [seq])

    def testGetPdbstr(self):
        cmd.pseudoatom()
        lines = cmd.get_pdbstr().splitlines()
        self.assertTrue(lines[0].startswith('HETATM'))
        self.assertTrue(lines[1].startswith('END'))

    def testGetSession(self):
        cmd.fragment('ala')
        x = cmd.count_atoms()
        s = cmd.get_session()
        cmd.reinitialize()
        cmd.set_session(s)
        self.assertEqual(['ala'], cmd.get_names())
        self.assertEqual(x, cmd.count_atoms())

    @testing.requires_version('2.5')
    def testGetSession25(self):
        cmd.set('pse_binary_dump', 0)
        cmd.set('pse_export_version', 0)
        s = cmd.get_session(binary=-1, version=-1)
        self.assertEqual(cmd.get_setting_int('pse_binary_dump'), 0)
        self.assertEqual(cmd.get_setting_float('pse_export_version'), 0.0)
        cmd.set_session(s)
        self.assertEqual(cmd.get_setting_int('pse_binary_dump'), 0)
        self.assertEqual(cmd.get_setting_float('pse_export_version'), 0.0)
        s = cmd.get_session(binary=True, version=1.2)
        self.assertEqual(cmd.get_setting_int('pse_binary_dump'), 0)
        self.assertEqual(cmd.get_setting_float('pse_export_version'), 0.0)
        cmd.set_session(s)
        self.assertEqual(cmd.get_setting_int('pse_binary_dump'), 1)
        self.assertEqual(cmd.get_setting_float('pse_export_version'), 1.2)

    @testing.requires_version('1.8.4')
    def testMultisave(self):
        names = ['ala', 'gly', 'his', 'arg']

        for name in names:
            cmd.fragment(name)

        sym1 = (1.2, 3.4, 5.6, 60.0, 70.0, 80.0)
        sym2 = (2.3, 4.5, 6.7, 90.0, 90.0, 90.0)
        sym3 = (3.4, 5.6, 7.8, 70.0, 70.0, 70.0)
        cmd.set_symmetry("gly", *sym1)
        cmd.set_symmetry("his", *sym2)
        cmd.set_symmetry("arg", *sym3)

        with testing.mktemp('.pdb') as filename:
            cmd.multisave(filename, names[0])                       # new
            cmd.multisave(filename, names[1])                       # new (overwrite)
            cmd.multisave(filename, ' '.join(names[2:]), append=1)  # append
            cmd.delete('*')
            cmd.load(filename)

        self.assertEqual(cmd.get_object_list(), names[1:])

        self.assertArrayEqual(cmd.get_symmetry("gly")[:6], sym1, delta=1e-4)
        self.assertArrayEqual(cmd.get_symmetry("his")[:6], sym2, delta=1e-4)
        self.assertArrayEqual(cmd.get_symmetry("arg")[:6], sym3, delta=1e-4)

    @testing.requires_version('2.1')
    def testMultifilesave(self):
        import glob

        for name in ['ala', 'gly', 'his', 'arg']:
            cmd.fragment(name)

        # multistate
        for i in range(2, 11):
            cmd.create('ala', 'ala', 1, i)

        for fmt in ['{}-{:02}.cif', '{}-{state}.cif']:
            with testing.mkdtemp() as dirname:
                cmd.multifilesave(os.path.join(dirname, fmt), 'ala', 0)
                filenames = [os.path.basename(p) for p in glob.glob(os.path.join(dirname, '*.cif'))]
                self.assertEqual(len(filenames), 10)
                self.assertTrue('ala-03.cif' in filenames)

        with testing.mkdtemp() as dirname:
            cmd.multifilesave(os.path.join(dirname, '{}.pdb'), 'a* g*')
            filenames_full = sorted(glob.glob(os.path.join(dirname, '*.pdb')))
            filenames = [os.path.basename(p) for p in filenames_full]
            self.assertEqual(filenames, ['ala.pdb', 'arg.pdb', 'gly.pdb'])

            cmd.delete('*')
            cmd.load(filenames_full[0])
            self.assertEqual(cmd.count_atoms(), 10)

            cmd.delete('*')
            cmd.load(filenames_full[1])
            self.assertEqual(cmd.count_atoms(), 24)

    @testing.foreach(
        ('0.5in', '0.25in', 200, (100, 50)),   # inch
        ('2.54cm', '1.27cm', 100, (100, 50)),  # centimeter
        (100, 0, -1, (100, 75)),               # px, only width
        (0, 75, -1, (100, 75)),                # px, only height
    )
    @testing.requires('no_edu') # ray
    @testing.requires_version('1.6.0')
    def testPng(self, w, h, dpi, size):
        with testing.mktemp('.png') as filename:
            cmd.png(filename, w, h, dpi, ray=1)
            self.assertEqual(size, Image.open(filename).size)

    @testing.requires_version('2.5')
    def testPngNoFile(self):
        self.ambientOnly()
        cmd.fragment('gly')
        cmd.show_as('spheres')
        cmd.color('yellow')
        cmd.zoom(complete=1)

        cmd.refresh()  # TODO should not be necessary, see also PYMOL-3328

        ncol, nrow = 120, 80
        buf = cmd.png(None, ncol, nrow)
        self.assertTrue(buf.startswith(b'\x89PNG'))
        import io
        img = self.get_imagearray(Image.open(io.BytesIO(buf)))
        self.assertEqual(img.shape[:2], (nrow, ncol))
        self.assertImageHasColor('yellow', img)

    # not supported in older versions: xyz (no ref)
    @testing.foreach('pdb', 'sdf', 'mol', 'mol2')
    def testSaveRef(self, format):
        # for rms_cur (not all formats save all identifiers)
        m = -1
        cmd.set('retain_order')

        cmd.fragment('ala', 'm1')
        cmd.copy('m2', 'm1')
        cmd.copy('m3', 'm1')

        cmd.rotate('y', 90, 'm2')
        cmd.align('m3', 'm2')

        # with ref=m3
        with testing.mktemp('.' + format) as filename:
            cmd.save(filename, 'm2', ref='m3')
            cmd.load(filename, 'm4')
            self.assertAlmostEqual(cmd.rms_cur('m4', 'm1', matchmaker=m), 0.00, delta=1e-2)
            self.assertAlmostEqual(cmd.rms_cur('m4', 'm2', matchmaker=m), 1.87, delta=1e-2)

        # without ref
        with testing.mktemp('.' + format) as filename:
            cmd.save(filename, 'm2')
            cmd.load(filename, 'm5')
            self.assertAlmostEqual(cmd.rms_cur('m5', 'm2', matchmaker=m), 0.00, delta=1e-2)
            self.assertAlmostEqual(cmd.rms_cur('m5', 'm1', matchmaker=m), 1.87, delta=1e-2)

    @testing.foreach(
            ('pdb',  1.2),
            ('sdf',  1.7),
            ('mol2', 1.7),
            ('xyz',  1.8),
            ('mol',  1.831),
            ('mae',  1.831),
            ('cif',  1.8),
    )
    def testSaveState(self, format, pymol_version):
        if pymol_version > testing.PYMOL_VERSION[1]:
            self.skipTest("version %f" % (pymol_version))

        # for rms_cur (not all formats save all identifiers)
        m = -1
        cmd.set('retain_order')

        # create a multistate object
        cmd.fragment('ala', 'm1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m1', 'm1', 1, 3)
        cmd.translate([5, 0, 0], 'm1', state=2)
        cmd.translate([0, 5, 0], 'm1', state=3)
        n_states = cmd.count_states('m1')

        with testing.mktemp('.' + format) as filename:
            if format == 'mae' and not pymol.invocation.options.incentive_product:
                format = 'cms'

            # explicit
            for state in range(1, n_states + 1):
                cmd.delete('m2')
                cmd.save(filename, 'm1', state=state)
                cmd.load(filename, 'm2', format=format)
                rms = cmd.rms_cur('m1', 'm2', state, 1, matchmaker=m)
                self.assertAlmostEqual(rms, 0.00, delta=1e-2)

            # current state
            for state in range(1, n_states + 1):
                cmd.frame(state)
                cmd.delete('m2')
                cmd.save(filename, 'm1')
                cmd.load(filename, 'm2', 1, format=format)
                rms = cmd.rms_cur('m1', 'm2', state, 1, matchmaker=m)
                self.assertAlmostEqual(rms, 0.00, delta=1e-2)

            if format in ('mol', 'cms'):
                # no multi support
                return

            # all states
            cmd.delete('m2')
            cmd.save(filename, 'm1', state=0)
            cmd.load(filename, 'm2', 1, discrete=1, multiplex=0)
            self.assertEqual(cmd.count_states('m2'), n_states)
            for state in range(1, n_states + 1):
                rms = cmd.rms_cur('m1', 'm2 and state %d' % state,
                        state, state, matchmaker=m)
                self.assertAlmostEqual(rms, 0.00, delta=1e-2)

    @testing.foreach(
            ('pdb',  1.2),
            ('sdf',  1.7),
            ('mol2', 1.7),
            ('xyz',  1.7),
            ('mol',  1.7),
            ('mae',  1.831),
            ('cif',  1.8),
    )
    def testSaveSelection(self, format, pymol_version):
        if pymol_version > testing.PYMOL_VERSION[1]:
            self.skipTest("version %f" % (pymol_version))

        cmd.fragment('trp', 'm1')
        cmd.fragment('glu', 'm2')

        n_O = cmd.count_atoms('elem O')
        n_N = cmd.count_atoms('elem N')

        with testing.mktemp('.' + format) as filename:
            if testing.PYMOL_VERSION_TUPLE < (1, 8):
                cmd.set('raise_exceptions', 0) # 1.7.6 save xyz doesn't set r=DEFAULT_SUCCESS
            cmd.save(filename, 'elem O+N')
            if testing.PYMOL_VERSION_TUPLE < (1, 8):
                cmd.set('raise_exceptions', 1)

            if format == 'mae' and not pymol.invocation.options.incentive_product:
                format = 'cms'

            cmd.delete('*')
            cmd.load(filename, 'm2', discrete=1, format=format) # avoid merging of atoms

        self.assertEqual(n_O, cmd.count_atoms('elem O'))
        self.assertEqual(n_N, cmd.count_atoms('elem N'))
        self.assertEqual(n_O + n_N, cmd.count_atoms())

    @testing.foreach('pdb', 'cif', 'mmtf')
    @testing.requires_version('2.5')
    def testSave_symmetry(self, format):
        sym1 = (1.2, 3.4, 5.6, 60.0, 70.0, 80.0, "P 3")
        cmd.fragment("gly", "m1")
        cmd.set_symmetry("m1", *sym1)

        with testing.mktemp('.' + format) as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.load(filename, "m2")

        sym = cmd.get_symmetry("m2")
        self.assertArrayEqual(sym[:6], sym1[:6], delta=1e-4)
        self.assertEqual(sym[6], sym1[6])

    # cmp_atom : compares all fields in Atom (see chempy/__init__.py)
    #            except the id (which is unique to the instance)
    def cmp_atom(self, selfobj,other):
        cmp = lambda a, b: a != b # py3k
        return \
                cmp(type(selfobj), type(other)) or \
                cmp(selfobj.segi, other.segi) or \
                cmp(selfobj.chain, other.chain) or \
                cmp(selfobj.resi_number, other.resi_number) or \
                cmp(selfobj.resi, other.resi) or \
                cmp(selfobj.resn, other.resn) or \
                cmp(selfobj.symbol, other.symbol) or \
                cmp(selfobj.name, other.name)

    # cmp_bond : compares all fields in Bond (see chempy/__init__.py)
    def cmp_bond(self, selfobj,other):
        cmp = lambda a, b: a != b # py3k
        return \
                cmp(selfobj.order, other.order) or \
                cmp(selfobj.index, other.index)

    def assertModelsAreSame(self, m1, m2):
        self.assertTrue(len(m1.atom) == len(m2.atom))
        idx = 0
        for m1atomidx in m1.atom:
            self.assertTrue(self.cmp_atom(m1atomidx, m2.atom[idx]) == 0)
            idx = idx + 1
        idx = 0
        for m1bondidx in m1.bond:
            self.assertTrue(self.cmp_bond(m1bondidx, m2.bond[idx]) == 0)
            idx = idx + 1

    @testing.requires_version('1.7.6')
    def testPSEBulkImport(self):
        cmd.load(self.datafile('1rx1_1766_bulk.pse.gz'))
        m1 = cmd.get_model()
        cmd.load(self.datafile('1rx1_176.pse.gz'))
        m2 = cmd.get_model()
        self.assertModelsAreSame(m1, m2)

    @testing.foreach.product((0, 1.7, 1.76, 1.8, 1.82, 1.9), (0, 1))
    @testing.requires_version('1.7.6.5')
    def testPSEBulkExportImport(self, pse_export_version, pse_binary_dump):
        with testing.mktemp('.pse') as filename:
            cmd.load(self.datafile("1oky-frag.pdb"))
            m1 = cmd.get_model()
            cmd.set("pse_export_version", pse_export_version)
            cmd.set("pse_binary_dump", pse_binary_dump)
            cmd.save(filename)
            cmd.reinitialize()
            cmd.load(filename)
            m2 = cmd.get_model()
            self.assertModelsAreSame(m1, m2)

    def testGetModelObjectName(self):
        cmd.load(self.datafile("1oky-frag.pdb"))
        cmd.load(self.datafile('1rna.cif'))

        m1 = cmd.get_model()
        cnt = cmd.count_atoms('%1oky-frag')

        self.assertEqual(m1.atom[0].model, '1oky-frag')
        self.assertEqual(m1.atom[-1].model, '1rna')
        self.assertEqual(m1.atom[cnt].model, '1rna')


    @testing.requires_version('2.1')
    def testMMTF(self):
        '''Styled MMTF export/import'''
        S = 0b10            # 1 << 1 spheres
        D = 0b1000000000    # 1 << 9 dots
        B = 2 # blue
        R = 4 # red

        cmd.fragment('gly')
        cmd.color(B)
        cmd.color(R, 'elem C')
        cmd.show_as('spheres')
        cmd.show_as('dots', 'elem C')

        with testing.mktemp('.mmtf') as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.load(filename)

        color_list = []
        reps_list = []

        cmd.iterate('*', 'color_list.append(color)', space=locals())
        cmd.iterate('*', 'reps_list.append(reps)', space=locals())

        self.assertEqual(color_list, [B, R, R, B, B, B, B])
        self.assertEqual(reps_list, [S, D, D, S, S, S, S])

    @testing.requires_version('2.1')
    def testMMTFExportSele(self):
        cmd.fab('ACDE')

        with testing.mktemp('.mmtf') as filename:
            cmd.save(filename, 'resn CYS+ASP')
            cmd.delete('*')
            cmd.load(filename)

        self.assertEqual(cmd.count_atoms(), 23)
        self.assertEqual(cmd.count_atoms('bound_to ASP/CG'), 3)
        self.assertEqual(cmd.get_model('ASP/CG ASP/OD1').bond[0].order, 2)
        self.assertEqual(cmd.get_model('ASP/CG ASP/OD2').bond[0].order, 1)
        self.assertEqual(cmd.get_model('CYS/C ASP/N').bond[0].order, 1)

    @testing.requires_version('2.1')
    def testMMTFExportEmpty(self):
        with testing.mktemp('.mmtf') as filename:
            cmd.save(filename)
            cmd.load(filename, 'm1')

        self.assertEqual(cmd.count_atoms(), 0)
        self.assertEqual(cmd.get_names(), ['m1'])

    @testing.requires_version('2.4')
    def testglTF(self):
        '''glTF export'''
        cmd.fragment('gly')

        with testing.mktemp('.gltf') as filename:
            self.assertEqual(cmd.save(filename), 0)
            cmd.delete('*')

    def testSaveAln(self):
        cmd.fab('ACDEFGH', 'm1')
        cmd.fab('ACDFGH', 'm2')
        cmd.align('m1', 'm2', cycles=0, object='aln')
        with testing.mktemp('.aln') as filename:
            cmd.save(filename)
            with open(filename) as handle:
                lines = list(handle)
        self.assertEqual(lines[0].split(), ['CLUSTAL'])
        self.assertEqual(lines[1].split(), [])
        self.assertEqual(lines[2].split(), ['m1', 'ACDEFGH'])
        self.assertEqual(lines[3].split(), ['m2', 'ACD-FGH'])

    @testing.requires_version('2.3')
    def testSaveAlnNucleic(self):
        cmd.load(self.datafile('1rna.cif'))
        cmd.create('m1', 'chain A and not resi 6-7')
        cmd.create('m2', 'chain B')
        cmd.alter('m1 and resi 1', 'resn = "DT"') # mimic DNA
        cmd.alter('m2 and resi 20', 'resn = "UNK"') # mimic nonstd residue
        cmd.align('m1', 'm2', cycles=0, object='aln')
        with testing.mktemp('.aln') as filename:
            cmd.save(filename)
            with open(filename) as handle:
                lines = list(handle)
        self.assertEqual(lines[0].split(), ['CLUSTAL'])
        self.assertEqual(lines[1].split(), [])
        self.assertEqual(lines[2].split(), ['m1', 'TUAUA--UAUAUAA'])
        self.assertEqual(lines[3].split(), ['m2', 'UUAUA?AUAUAUAA'])

    @testing.requires_version('2.0')
    def testAssingAtomTypes(self):
        cmd.fragment('his')
        cmd.alter('all', 'text_type = "none"')
        mytypes = {}
        cmd.iterate('all', 'mytypes[name] = text_type', space=locals())
        self.assertEqual(mytypes['NE2'], 'none')
        pymol.exporting.assign_atom_types('all', 'mol2')
        cmd.iterate('all', 'mytypes[name] = text_type', space=locals())
        self.assertEqual(mytypes['NE2'], 'N.pl3')
        self.assertEqual(mytypes['N'], 'N.3')
        self.assertEqual(mytypes['CG'], 'C.2')
