'''
Testing: pymol.editing
'''

from pymol import cmd, testing, stored
from pymol import CmdException

def get_coord_list(selection, state=1):
    return cmd.get_model(selection, state).get_coord_list()

def get_atom_names(selection='all'):
    return [a.name for a in cmd.get_model(selection).atom]

class TestEditing(testing.PyMOLTestCase):

    alter_names_atomic = {
        str: ('segi', 'chain', 'resn', 'resi', 'name', 'alt', 'elem', 'text_type'),
        int: ('formal_charge', 'numeric_type', 'ID'),
        float: ('q', 'b', 'partial_charge', 'vdw'),
    }

    alter_names_atomic_special = {
        str: ('model', 'type'),
        int: ('resv',),
    }

    alter_names_state = {
        float: ('x', 'y', 'z'),
    }

    def testIterate(self):
        cmd.pseudoatom('m1')

        for ptype, names in (
                list(self.alter_names_atomic.items()) +
                list(self.alter_names_atomic_special.items())):
            stored.v_names = None
            cmd.iterate('all', 'stored.v_names = [' + ','.join(names) + ']')

            self.assertNotEqual(stored.v_names, None)
            for v in stored.v_names:
                self.assertTrue(isinstance(v, ptype))

    @testing.foreach(cmd.iterate, cmd.alter)
    def testIterateLocalVars(self, func):
        cmd.pseudoatom('m1')
        stored.v = None
        func('all', 'foobar = model;stored.v = foobar')
        self.assertEqual(stored.v, 'm1')

    @testing.requires_version("2.5")
    def testIterateCallback(self):
        cmd.fragment('gly')
        atoms = []
        cmd.iterate('not hydro', atoms.append)
        self.assertEqual([a.name for a in atoms], ['N', 'CA', 'C', 'O'])
        self.assertEqual([a.protons for a in atoms], [7, 6, 6, 8])

    @testing.requires_version("2.6")
    def testIterate__explicit_degree(self):
        cmd.fragment("ala")
        result = []
        cmd.iterate("all", "result.append(explicit_degree)", space=locals())
        self.assertEqual(result, [2, 4, 2, 1, 4, 1, 1, 1, 1, 1])

    @testing.requires_version("2.6")
    def testIterate__explicit_valence(self):
        cmd.fragment("ala")
        result = []
        cmd.iterate("all", "result.append(explicit_valence)", space=locals())
        self.assertEqual(result, [2, 4, 3, 2, 4, 1, 1, 1, 1, 1])

    @testing.requires_version("2.6")
    def testIterate__explicit_valence__aromatic(self):
        cmd.fragment("benzene")
        cmd.valence("aromatic", "elem C", "elem C")
        result = []
        cmd.iterate("all", "result.append(explicit_valence)", space=locals())
        self.assertEqual(result, [4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1, 1])

    def testIterateState(self):
        cmd.fragment('ala')
        cmd.create('ala', 'ala', 1, 2)
        self.assertEqual(2, cmd.count_states())

        v_count = cmd.count_atoms('all')
        expr = 'v_xyz.append(((model,index), (x,y,z)))'

        # current state
        v_xyz = []
        cmd.iterate_state(-1, 'all', expr, space=locals())
        self.assertEqual(len(v_xyz), v_count)

        # all states
        v_xyz = []
        cmd.iterate_state(0, 'all', expr, space=locals())
        self.assertEqual(len(v_xyz), v_count * 2)

        # atomic=0
        stored.v_xyz = []
        cmd.iterate_state(1, 'all', 'stored.v_xyz.append((x,y,z))', atomic=0)
        self.assertEqual(len(stored.v_xyz), v_count)

        space = {'self': self, 'NameError': NameError, 'v_list': []}
        cmd.iterate_state(1, 'all',
                'v_list.append(self.assertRaises(NameError, lambda: (model, index)))',
                space=space)
        self.assertEqual(len(space['v_list']), v_count)

    def testAlter(self):
        cmd.pseudoatom('m1')

        for ptype, names in self.alter_names_atomic.items():
            stored.v_mock = tuple(map(ptype, range(len(names))))
            names_fmt = '(' + ','.join(names) + ',)'

            if not (cmd.alter('all', names_fmt + ' = stored.v_mock') and
                    cmd.iterate('all', 'stored.v_names = ' + names_fmt)):
                raise UserWarning('iterate failed')

            self.assertEqual(stored.v_names, stored.v_mock)

    def testAlterState(self):
        cmd.fragment('ala')
        cmd.create('ala', 'ala', 1, 2)
        v_count = cmd.count_atoms('all')

        v_mock = [list(map(float, list(range(i*3, (i+1)*3)))) for i in range(v_count)]

        v_xyz_1_pre = cmd.get_model('all', state=1).get_coord_list()
        v_xyz_2_pre = cmd.get_model('all', state=2).get_coord_list()

        cmd.alter_state(2, 'all', '(x,y,z) = next(xyz_iter)',
                space={'xyz_iter': iter(v_mock), 'next': next})

        v_xyz_1_post = cmd.get_model('all', state=1).get_coord_list()
        v_xyz_2_post = cmd.get_model('all', state=2).get_coord_list()

        self.assertEqual(v_xyz_1_post, v_xyz_1_pre)
        self.assertEqual(v_xyz_2_post, v_mock)

    def test_alter_list(self):
        cmd.fragment('gly')
        cmd.alter_list('gly', [[i+1, 'name = "X%d"' % i] for i in range(7)])
        name_list = []
        cmd.iterate('gly', 'name_list.append(name)', space=locals())
        self.assertEqual(name_list, ['X%d' % i for i in range(7)])

    @testing.requires_version('2.5')
    def test_alter_exceptions(self):
        cmd.fragment('gly')
        with self.assertRaisesRegex(Exception,
                                    'Use alter/alter_state to modify values'):
            cmd.iterate('all', 'b = "abc"')
        with self.assertRaisesRegex(ValueError, 'float'):
            cmd.alter('all', 'b = "abc"')
        with self.assertRaises(IndexError):
            cmd.iterate('all', 'name[100]')

    def test_attach(self):
        cmd.pseudoatom()
        cmd.edit('first all')
        cmd.attach('C', 1, 1)
        self.assertEqual(2, cmd.count_atoms())
        cmd.attach('C', 1, 1)
        self.assertEqual(3, cmd.count_atoms())
        self.assertEqual(['C01', 'C02', 'PS1'], sorted(get_atom_names()))

    @testing.requires_version('2.4')
    def test_add_bond(self):
        cmd.pseudoatom('m1', pos=(0,0,0))
        cmd.pseudoatom('m1', pos=(1,0,0))
        cmd.pseudoatom('m1', pos=(1,1,0))
        cmd.add_bond('m1', 1, 2) # 1-indexed
        cmd.add_bond('m1', 2, 3, order=2)
        bonds = cmd.get_bonds('m1') # 0-indexed
        self.assertEqual(bonds, [(0, 1, 1), (1, 2, 2)])

    @testing.requires_version('2.4')
    def test_rebond(self):
        cmd.pseudoatom('m1', pos=(0, 0, 0), vdw=0.8)
        cmd.pseudoatom('m1', pos=(1, 0, 0), vdw=0.8)
        cmd.pseudoatom('m1', pos=(1, 1, 0), vdw=0.8)
        cmd.add_bond('m1', 1, 3) # this bond will be removed by "rebond"
        cmd.rebond('m1')
        bonds = cmd.get_bonds('m1') # 0-indexed
        self.assertEqual(bonds, [(0, 1, 1), (1, 2, 1)])

    def test_bond(self):
        cmd.pseudoatom('m1', pos=(0,0,0))
        cmd.pseudoatom('m1', pos=(1,0,0))
        cmd.pseudoatom('m1', pos=(1,1,0))

        cmd.bond('m1`1', 'm1`2')
        count = cmd.count_atoms('(m1`1) extend 1')
        self.assertEqual(count, 2)

        cmd.unbond('m1`1', 'm1`2')
        count = cmd.count_atoms('(m1`1) extend 1')
        self.assertEqual(count, 1)

    def test_cycle_valence(self):
        cmd.fragment('gly')
        cmd.edit('ID 0', 'ID 1')

        cmd.cycle_valence()
        self.assertEqual(4, cmd.get_model('pkbond').bond[0].order)

        cmd.cycle_valence()
        self.assertEqual(2, cmd.get_model('pkbond').bond[0].order)

        cmd.cycle_valence()
        self.assertEqual(3, cmd.get_model('pkbond').bond[0].order)

    def test_deprotect(self):
        # see test_protect
        pass

    def test_drag(self):
        cmd.drag
        self.skipTest("TODO")

    def test_dss(self):
        ss_list = []
        cmd.fab('A' * 6, ss=1)
        cmd.dss()
        cmd.iterate('2-5/CA', 'ss_list.append(ss)', space=locals())
        self.assertEqual(ss_list, ['H', 'H', 'H', 'H'])

    def test_edit(self):
        cmd.fragment('gly')
        cmd.edit('ID 0', 'ID 1', 'ID 2','ID 3')
        names = cmd.get_names('public_selections')
        self.assertEqual(names, ['pk1', 'pk2', 'pk3', 'pk4', 'pkset', 'pkmol'])

    def test_fix_chemistry(self):

        # Baseline bond order
        cmd.fragment('arg')
        # list(set()) to remove possibly duplicate bond info
        bonds = list(set(cmd.get_bonds('arg')))
        total_order = sum(t[2] for t in bonds)
        self.assertEqual(total_order, 25)

        # Change bond orders
        cmd.bond('arg`9', 'arg`11', 2)
        bonds = list(set(cmd.get_bonds('arg')))
        total_order = sum(t[2] for t in bonds)
        self.assertEqual(total_order, 27)

        # Restore with fix_chemistry
        cmd.fix_chemistry('arg')
        bonds = list(set(cmd.get_bonds('arg')))
        total_order = sum(t[2] for t in bonds)
        self.assertEqual(total_order, 25)


    def test_flag(self):
        cmd.fab("ACE")
        natoms = cmd.count_atoms()
        natoms_cys = cmd.count_atoms('resn CYS')
        assert natoms == 36
        assert natoms_cys == 11
        cmd.alter("all", "flags = 0")
        for flag in range(32):
            sele = "flag {}".format(flag)
            self.assertEqual(cmd.count_atoms(sele), 0)
            cmd.flag(flag, "index 1-{}".format(flag + 1), "set")
            self.assertEqual(cmd.count_atoms(sele), flag + 1)
            cmd.flag(flag, "index {}-{}".format(flag + 1, flag + 2), "reset")
            self.assertEqual(cmd.count_atoms(sele), 2)
            cmd.flag(flag, "index 1-{}".format(flag + 1), "clear")
            self.assertEqual(cmd.count_atoms(sele), 1)
            cmd.flag(flag, "last all", "set")
            self.assertEqual(cmd.count_atoms(sele), 2)
        cmd.alter("all", "flags = 0")
        for flag, flagname in [
            (0, 'focus'),
            (1, 'free'),
            (2, 'restrain'),
            (3, 'fix'),
            (4, 'exclude'),
            (5, 'study'),
            (24, 'exfoliate'),  # deprecated
            (25, 'ignore'),
            (26, 'no_smooth'),
        ]:
            sele = "flag {}".format(flag)
            self.assertEqual(cmd.count_atoms(sele), 0)
            cmd.flag(flagname, 'resn CYS')
            self.assertEqual(cmd.count_atoms(sele), natoms_cys)

        cmd.alter("all", "flags = 1 << 25")
        self.assertEqual(cmd.count_atoms("flag 25"), natoms)
        # number as string
        cmd.flag("25", "index 11-15")  # reset
        self.assertEqual(cmd.count_atoms("flag 25"), 5)
        # abbreviated action
        cmd.flag("25", "index 15-20", "s")  # set
        self.assertEqual(cmd.count_atoms("flag 25"), 10)

    @testing.requires_version('2.5')
    def test_flag_error_handling(self):
        with self.assertRaisesRegex(CmdException, "out of range"):
            cmd.flag(-1, "all")
        with self.assertRaisesRegex(CmdException, "out of range"):
            cmd.flag(32, "all")
        with self.assertRaisesRegex(cmd.QuietException, "unknown action"):
            cmd.flag(1, "all", "foo")
        with self.assertRaisesRegex(CmdException, "Invalid selection name"):
            cmd.flag(1, "foo")

    # 2.1 for get_bonds
    @testing.requires_version('2.1')
    def test_fuse(self):
        cmd.fragment('ala')
        cmd.fragment('gly')
        # non-hydrogens (N has an open valence, O actually not)
        cmd.fuse("gly and elem N", "ala and elem O")
        self.assertEqual(17, cmd.count_atoms('ala'))
        # hydrogens
        cmd.fuse("/gly///GLY/3HA", "/ala///ALA/3HB")
        self.assertEqual(22, cmd.count_atoms("ala"))
        # mode 3: don't create a bond
        cmd.fuse("gly", "ala", mode=3)
        self.assertEqual(29, cmd.count_atoms("ala"))
        # move=0
        cmd.fuse("/gly///GLY/C", "/ala///GLY/C03", move=0)
        self.assertEqual(36, cmd.count_atoms("ala"))
        self.assertEqual(cmd.get_bonds("ala"), [
            (0, 1, 1),
            (0, 5, 1),
            (1, 2, 1),
            (1, 4, 1),
            (1, 6, 1),
            (2, 3, 2),
            (3, 9, 1),
            (4, 13, 1),
            (4, 7, 1),
            (4, 8, 1),
            (9, 16, 1),
            (9, 25, 1),
            (10, 13, 1),
            (10, 26, 1),
            (11, 14, 1),
            (11, 27, 1),
            (12, 15, 1),
            (12, 28, 1),
            (13, 18, 1),
            (13, 29, 1),
            (14, 19, 1),
            (14, 30, 1),
            (14, 31, 1),
            (15, 20, 1),
            (15, 32, 1),
            (15, 33, 1),
            (16, 17, 1),
            (16, 34, 1),
            (16, 35, 1),
            (17, 21, 2),
            (18, 20, 1),
            (18, 22, 2),
            (19, 23, 2),
            (20, 24, 2),
        ])
        self.assertArrayEqual(cmd.get_coordset("ala"), [
            [-0.68, -1.23, -0.49],
            [-0.0, 0.06, -0.49],
            [-0.51, 0.86, 0.73],
            [1.5, -0.11, -0.49],
            [2.03, -1.23, -0.5],
            [-1.6, 1.01, 0.69],
            [-0.28, 0.34, 1.68],
            [-0.13, -2.16, -0.49],
            [-0.27, 0.6, -1.42],
            [1.24, -2.38, -0.51],
            [-0.1, -1.95, -0.91],
            [-1.12, -2.31, 0.15],
            [-0.8, -2.9, 1.19],
            [1.42, -2.89, 0.42],
            [-0.35, -2.43, -1.81],
            [-0.11, -0.85, -1.02],
            [1.55, 2.1, 1.1],
            [0.13, 2.26, 0.81],
            [-0.57, 3.06, 1.88],
            [0.04, 3.5, 2.86],
            [2.0, 2.52, 1.99],
            [0.02, 2.8, -0.15],
            [-1.19, 0.2, -0.21],
            [0.23, 0.32, -0.5],
            [1.06, -0.39, 0.54],
            [0.55, -0.97, 1.5],
            [-1.56, -0.33, 0.66],
            [0.48, 1.34, -0.51],
            [0.43, -0.16, -1.48],
            [-1.19, 0.2, -0.21],
            [0.23, 0.32, -0.5],
            [1.06, -0.39, 0.54],
            [0.55, -0.97, 1.5],
            [-1.56, -0.33, 0.66],
            [0.48, 1.34, -0.51],
            [0.43, -0.16, -1.48],
        ], 1e-2)
        names = []
        cmd.iterate("ala", "names.append(name)", space=locals())
        self.assertEqual(names, [
            'N', 'CA', 'C', 'O', 'CB', 'H', 'HA', '1HB', '2HB', 'N', 'N01',
            'N02', 'N03', 'C02', 'C04', 'C06', 'CA', 'C', 'C03', 'C05', 'C07',
            'O', 'O04', 'O06', 'O08', 'H', 'H05', 'H08', 'H11', 'H07', 'H09',
            'H10', 'H12', 'H13', '3HA', 'HA'
        ])

    def test_get_editor_scheme(self):
        cmd.get_editor_scheme
        self.skipTest("TODO")

    def test_h_add(self):
        cmd.fragment('gly')
        cmd.h_add()
        self.assertEqual(5, cmd.count_atoms('hydro'))

    @testing.requires_version('2.1')
    def test_h_add_state(self):
        nheavy = 4
        nhydro = 5
        nfull = nheavy + nhydro

        cmd.fragment('gly', 'm1')
        cmd.remove('hydro')
        self.assertEqual(nheavy, cmd.count_atoms())
        cmd.h_add()
        self.assertEqual(nfull, cmd.count_atoms())

        # multi-state
        cmd.remove('hydro')
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m1', 'm1', 1, 3)
        cmd.h_add(state=2)
        self.assertEqual(nfull, cmd.count_atoms())
        self.assertEqual(nheavy, cmd.count_atoms('state 1'))
        self.assertEqual(nfull, cmd.count_atoms('state 2'))
        self.assertEqual(nheavy, cmd.count_atoms('state 3'))

        # discrete multi-state
        cmd.remove('hydro')
        cmd.create('m2', 'm1', 1, 1, discrete=1)
        cmd.create('m2', 'm2', 1, 2, discrete=1)
        cmd.create('m2', 'm2', 1, 3, discrete=1)
        self.assertEqual(nheavy * 3, cmd.count_atoms('m2'))
        cmd.h_add('m2 & state 2') # TODO , state=2)
        self.assertEqual(nfull + nheavy * 2, cmd.count_atoms('m2'))
        self.assertEqual(nheavy, cmd.count_atoms('m2 & state 1'))
        self.assertEqual(nfull, cmd.count_atoms('m2 & state 2'))
        self.assertEqual(nheavy, cmd.count_atoms('m2 & state 3'))
        cmd.h_add('m2')
        self.assertEqual(nfull * 3, cmd.count_atoms('m2'))


    def test_h_fill(self):
        cmd.fragment('gly')
        cmd.edit('elem N')
        cmd.h_fill()
        self.assertEqual(4, cmd.count_atoms('hydro'))

    def test_h_fix(self):
        cmd.h_fix
        self.skipTest("TODO")

    def test_invert(self):
        cmd.fragment('ala')
        xyzfix = get_coord_list('ID 1-7')
        xyzmov = get_coord_list('ID 0+8+9')
        cmd.edit('ID 1', 'ID 2', 'ID 3')
        cmd.invert()
        self.assertEqual(xyzfix, get_coord_list('ID 1-7'))
        self.assertNotEqual(xyzmov, get_coord_list('ID 0+8+9'))

    def test_map_double(self):
        cmd.map_double
        self.skipTest("TODO")

    def test_map_halve(self):
        cmd.map_halve
        self.skipTest("TODO")

    def test_map_set(self):
        cmd.map_set
        self.skipTest("TODO")

    def test_map_set_border(self):
        cmd.map_set_border
        self.skipTest("TODO")

    def test_map_trim(self):
        cmd.map_trim
        self.skipTest("TODO")

    def test_matrix_copy(self):
        cmd.fragment('ala')
        cmd.fragment('gly')
        cmd.rotate('x', 90, 'none', object='ala')
        cmd.matrix_copy('ala', 'gly')
        self.assertArrayEqual(cmd.get_object_matrix('gly'),
                (1.0, 0.0, 0.0, 0.0,
                 0.0, 0.0,-1.0, 0.0,
                 0.0, 1.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 1.0), 1e-4)
        cmd.matrix_reset('gly', mode=1)
        self.assertArrayEqual(cmd.get_object_matrix('gly'),
                (1.0, 0.0, 0.0, 0.0,
                 0.0, 1.0, 0.0, 0.0,
                 0.0, 0.0, 1.0, 0.0,
                 0.0, 0.0, 0.0, 1.0), 1e-4)

    def test_matrix_reset(self):
        # see test_matrix_copy
        pass

    def test_protect(self):
        cmd.pseudoatom('m1', pos=[0.,0.,0.])
        cmd.pseudoatom('m1', pos=[1.,0.,0.])

        cmd.protect('m1`1')
        cmd.translate([0.,1.,0.])
        self.assertEqual([0.,0.,0.], cmd.get_atom_coords('m1`1'))
        self.assertEqual([1.,1.,0.], cmd.get_atom_coords('m1`2'))

        cmd.deprotect()
        cmd.translate([0.,0.,1.])
        self.assertEqual([0.,0.,1.], cmd.get_atom_coords('m1`1'))
        self.assertEqual([1.,1.,1.], cmd.get_atom_coords('m1`2'))

    def test_push_undo(self):
        cmd.push_undo
        self.skipTest("TODO")

    def test_redo(self):
        cmd.redo
        self.skipTest("TODO")

    def test_reference(self):
        # undocumented
        cmd.reference
        self.skipTest("TODO")

    def test_remove(self):
        cmd.fragment('ala')
        cmd.remove('elem C')
        self.assertEqual(7, cmd.count_atoms())

    @testing.foreach((0, 7), (1, 3))
    def test_remove_picked(self, hydro, count):
        cmd.fragment('ala')
        cmd.edit('ID 1', 'ID 2', 'ID 3')
        cmd.remove_picked(hydro)
        self.assertEqual(count, cmd.count_atoms())

    def test_rename(self):
        cmd.fragment('ala')
        n_atoms = cmd.count_atoms()
        v = cmd.alter('elem C', 'name="C"')
        self.assertEqual(v, 3)
        self.assertEqual(['C', 'C', 'C'], get_atom_names('elem C'))
        cmd.rename('ala')
        self.assertEqual(['C', 'C01', 'C02'], sorted(get_atom_names('elem C')))
        cmd.alter('index 1-3', 'name="X"')
        cmd.alter('index 4', 'name="Y"')
        self.assertEqual(3, cmd.count_atoms('name X'))
        self.assertEqual(1, cmd.count_atoms('name Y'))
        cmd.rename('ala', force=1)
        atom_names = get_atom_names()
        self.assertTrue('X' not in atom_names)
        self.assertTrue('Y' not in atom_names)
        self.assertEqual(n_atoms, len(set(atom_names)))

    def test_replace(self):
        # NOTES:
        # - doc says "Immature functionality"
        # - will not preserve ID
        # - will not re-add hydrogens with h_fill=1
        cmd.fragment('ala')
        cmd.edit('elem O')
        cmd.replace('S', 0, 0, 0)
        self.assertEqual(0, cmd.count_atoms('elem O'))
        self.assertEqual(1, cmd.count_atoms('elem S'))
        self.assertEqual(1, cmd.count_atoms('name S01'))
        self.assertEqual(0, cmd.count_atoms('name S02'))
        cmd.edit('elem N')
        cmd.replace('S', 0, 0, 0)
        self.assertEqual(2, cmd.count_atoms('elem S'))
        self.assertEqual(1, cmd.count_atoms('name S02'))

    def test_rotate(self):
        from chempy import cpv
        from math import pi

        get_coord_list = lambda s: cmd.get_model(s).get_coord_list()

        cmd.fragment('ala')
        cmd.select('s1', 'hydro')

        c1_hydro = get_coord_list('s1')
        c1_other = get_coord_list('not s1')
        c1_hydro_rot = cpv.transform_array(
                cpv.rotation_matrix(pi * 0.25, [0,0,1]), c1_hydro)

        cmd.rotate('z', 45, 's1', camera=0, origin=(0,0,0))

        c2_hydro = get_coord_list('s1')
        c2_other = get_coord_list('not s1')

        self.assertArrayEqual(c2_other, c1_other)
        self.assertArrayEqual(c2_hydro, c1_hydro_rot, 0.001)

    def test_sculpt_activate(self):
        cmd.sculpt_activate
        self.skipTest("TODO")

    def test_sculpt_deactivate(self):
        cmd.sculpt_deactivate
        self.skipTest("TODO")

    def test_sculpt_iterate(self):
        cmd.sculpt_iterate
        self.skipTest("TODO")

    def test_sculpt_purge(self):
        cmd.sculpt_purge
        self.skipTest("TODO")

    def test_set_dihedral(self):
        cmd.fragment('ala')
        angle = 45.0
        atoms = ('ID 7', 'ID 2', 'ID 1', 'ID 9')
        cmd.set_dihedral(*atoms, angle=angle)
        v = cmd.get_dihedral(*atoms)
        self.assertAlmostEqual(v, angle, 4)

    def test_set_geometry(self):
        cmd.set_geometry
        self.skipTest("TODO")

    def test_set_name(self):
        cmd.pseudoatom('m1')
        cmd.set_name('m1', 'm2')
        self.assertEqual(['m2'], cmd.get_names())

    def test_set_object_color(self):
        self.skipTest("TODO")

        cmd.fragment('gly')
        cmd.set_object_color('gly', 'blue')
        self.assertEqual('blue', cmd.get('color', 'gly'))

    def test_set_object_ttt(self):
        M = [1.0, 0.0, 0.0, 0.0,
             0.0, 0.0,-1.0, 0.0,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 1.0]
        cmd.pseudoatom('m1')
        cmd.set_object_ttt('m1', M)
        self.assertArrayEqual(M, cmd.get_object_matrix('m1'), 1e-4)

    def test_set_symmetry(self):
        sym = [68.7, 126.8, 184.0, 90.0, 90.0, 90.0, 'P 21 21 21']
        cmd.pseudoatom('m1')
        cmd.set_symmetry('m1', *sym)
        v = cmd.get_symmetry('m1')
        self.assertEqual(v[-1], sym[-1])
        self.assertArrayEqual(v[:-1], sym[:-1], 1e-4)

        cmd.pseudoatom('m2')
        cmd.symmetry_copy('m1', 'm2')
        v = cmd.get_symmetry('m2')
        self.assertEqual(v[-1], sym[-1])
        self.assertArrayEqual(v[:-1], sym[:-1], 1e-4)

    @testing.requires_version('2.5')
    def test_set_symmetry__state(self):
        sym1 = [68.7, 126.8, 184.0, 90.0, 90.0, 90.0, 'P 21 21 21']
        sym2 = [12, 34, 56, 90.0, 90.0, 90.0, 'P 1']
        cmd.pseudoatom('m1')
        cmd.create('m1', 'm1', 1, 2)
        assert cmd.count_states('m1') == 2
        cmd.set_symmetry('m1', *sym1, state=1)
        cmd.set_symmetry('m1', *sym2, state=2)
        for globalstate, state, sym in [
            (2, 1, sym1),
            (1, 2, sym2),
            (1, -1, sym1),  # current state
            (2, -1, sym2),  # current state
        ]:
            cmd.frame(globalstate)
            v = cmd.get_symmetry('m1', state=state)
            self.assertEqual(v[-1], sym[-1])
            self.assertArrayEqual(v[:-1], sym[:-1], 1e-4)

        cmd.pseudoatom('m2', state=1)
        cmd.create('m2', 'm2', 1, 2)
        assert cmd.count_states('m2') == 2
        cmd.symmetry_copy('m1', 'm2', 1, 2)
        cmd.symmetry_copy('m1', 'm2', 2, 1)
        self.assertEqual(cmd.get_symmetry('m2', 1)[-1], sym2[-1])
        self.assertEqual(cmd.get_symmetry('m2', 2)[-1], sym1[-1])

        # rendering
        self.ambientOnly()
        cmd.disable('m2')
        cmd.color('blue', 'm1')
        cmd.set('cgo_line_width', 10)
        cmd.show_as('cell', 'm1')
        cmd.frame(1)
        cmd.set_view ((
            0.864273727,    0.052246489,   -0.500302136,
            0.200916424,    0.875955939,    0.438560009,
            0.461155534,   -0.479554355,    0.746567726,
            0.000069976,    0.000004172, -809.162963867,
            71.266143799,   17.057806015,   37.736812592,
            560.533935547, 1097.791748047,  -20.0))
        img1 = self.get_imagearray(width=100, height=100)
        self.assertImageHasColor('blue', img1)
        cmd.frame(2)
        img2 = self.get_imagearray(width=100, height=100)
        self.assertImageHasColor('blue', img2)
        self.assertImageNotEqual(img1, img2)

    @testing.requires_version('1.7.3.0')
    def test_set_state_order(self):
        import numpy

        cmd.fragment('ala', 'm1')
        for i in (2, 3):
            cmd.create('m1', 'm1', i - 1, i);
            cmd.rotate('x', 10.0, 'm1', i)

        coords1 = cmd.get_coordset('m1', 1)
        coords3 = cmd.get_coordset('m1', 3)

        self.assertEqual(coords1.shape, (10, 3))
        self.assertFalse(numpy.allclose(coords1, coords3))

        cmd.set_state_order('m1', (3, 2, 1))

        self.assertTrue(numpy.allclose(coords1, cmd.get_coordset('m1', 3)))
        self.assertTrue(numpy.allclose(coords3, cmd.get_coordset('m1', 1)))

    def test_set_title(self):
        text = 'foo'
        cmd.pseudoatom('m1')
        cmd.set_title('m1', 1, text)
        self.assertEqual(cmd.get_title('m1', 1), text)

    @testing.requires_version('2.4')
    def test_set_title_current_state(self):
        cmd.pseudoatom('m1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m1', 'm1', 1, 3)
        cmd.create('m1', 'm1', 1, 4)
        # global state
        cmd.frame(2)
        cmd.set_title('m1', -1, "foo")
        self.assertEqual(cmd.get_title('m1', -1), "foo") # current state
        self.assertEqual(cmd.get_title('m1', -2), "foo") # effective state (same as current)
        self.assertEqual(cmd.get_title('m1', 1), "")
        self.assertEqual(cmd.get_title('m1', 2), "foo")
        self.assertEqual(cmd.get_title('m1', 3), "")
        self.assertEqual(cmd.get_title('m1', 4), "")
        # object state
        cmd.set('state', 3, 'm1')
        cmd.set_title('m1', -1, "bar")
        self.assertEqual(cmd.get_title('m1', -1), "bar") # current state
        self.assertEqual(cmd.get_title('m1', -2), "bar") # effective state (same as current)
        self.assertEqual(cmd.get_title('m1', 1), "")
        self.assertEqual(cmd.get_title('m1', 2), "foo")
        self.assertEqual(cmd.get_title('m1', 3), "bar")
        self.assertEqual(cmd.get_title('m1', 4), "")

    def test_smooth(self):
        pdbfile = self.datafile("sampletrajectory.pdb")
        dcdfile = self.datafile("sampletrajectory.dcd")

        cmd.load(pdbfile)
        cmd.load_traj(dcdfile, state=0)

        cmd.select('ref_atom','/sampletrajectory///LYS`244/NZ')

        pre_1_xyz = []
        pre_5_xyz = []
        post_1_xyz = []
        post_5_xyz = []

        cmd.iterate_state(1, 'ref_atom', 'pre_1_xyz.append([x,y,z])', space=locals(), atomic=0)
        cmd.iterate_state(5, 'ref_atom', 'pre_5_xyz.append([x,y,z])', space=locals(), atomic=0)

        cmd.smooth('sampletrajectory',10,10)

        cmd.iterate_state(1, 'ref_atom', 'post_1_xyz.append([x,y,z])', space=locals(), atomic=0)
        cmd.iterate_state(5, 'ref_atom', 'post_5_xyz.append([x,y,z])', space=locals(), atomic=0)

        self.assertArrayEqual(pre_1_xyz[0], post_1_xyz[0], delta=0.01)
        self.assertArrayNotEqual(pre_5_xyz[0], post_5_xyz[0], delta=0.01)

        ref_values_5 = [43.79, 50.05, 14.93]
        self.assertArrayEqual(ref_values_5, post_5_xyz[0], delta=0.01)

        self.assertEqual(cmd.count_states('sampletrajectory'),11)

    def test_sort(self):
        cmd.pseudoatom('m1', name='PS2')
        cmd.pseudoatom('m1', name='PS1')
        cmd.alter('all', 'name = name_list.pop()', space={'name_list': ['PS1', 'PS2']})
        v = [a.name for a in cmd.get_model().atom]
        self.assertEqual(['PS2', 'PS1'], v)
        cmd.sort()
        v = [a.name for a in cmd.get_model().atom]
        self.assertEqual(['PS1', 'PS2'], v)

    def test_split_states(self):
        cmd.fragment('ala', 'm1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m1', 'm1', 1, 3)
        cmd.split_states('m1', 1, 2)
        self.assertItemsEqual(['m1', 'm1_0001', 'm1_0002'], cmd.get_names())

    def test_symmetry_copy(self):
        # see test_set_symmetry
        pass

    def test_torsion(self):
        cmd.fragment('ala')
        delta = 10
        atoms = ('ID 7', 'ID 2', 'ID 1', 'ID 9')
        d1 = cmd.get_dihedral(*atoms)
        cmd.edit(*atoms[1:3])
        cmd.torsion(delta)
        d2 = cmd.get_dihedral(*atoms)
        self.assertAlmostEqual(d1 + delta, d2, 4)

    def test_transform_object(self):
        cmd.transform_object
        self.skipTest("TODO")

    def test_transform_selection(self):
        cmd.transform_selection
        self.skipTest("TODO")

    def test_translate(self):
        # see test_protect
        pass

    def test_translate_atom(self):
        cmd.translate_atom
        self.skipTest("TODO")

    def test_unbond(self):
        # see test_bond
        pass

    def test_unpick(self):
        cmd.pseudoatom('m1')
        cmd.edit('m1')
        cmd.unpick()
        self.assertTrue('pk1' not in cmd.get_names('selections'))

    def test_update(self):
        # 3 states
        cmd.fragment('gly', 'm1')
        cmd.create('m1', 'm1', 1, 2)
        cmd.create('m1', 'm1', 1, 3)

        # second object, 90 degree rotates
        cmd.copy('m2', 'm1')
        cmd.rotate('x', 90, '(m2)', state=0)

        # reference coordsets
        cs = cmd.get_coordset
        cs1 = cs('m1', 1)
        cs2 = cs('m2', 1)

        # m2/3 will change (pre-check)
        self.assertArrayEqual(cs2, cs('m2', 3))
        self.assertArrayNotEqual(cs1, cs('m2', 3))

        # update explicit state
        cmd.update('m2', 'm1', 3, 2)

        # m2/3 has changed
        self.assertArrayEqual(cs1, cs('m2', 3))
        self.assertArrayNotEqual(cs2, cs('m2', 3))

        # these haven't changed
        self.assertArrayEqual(cs2, cs('m2', 1))
        self.assertArrayEqual(cs2, cs('m2', 2))

        # reset m2/3
        cmd.load_coordset(cs2, 'm2', 3)
        self.assertArrayEqual(cs2, cs('m2', 3))

        # update all states
        cmd.update('m2', 'm1', 0, 0)
        self.assertArrayEqual(cs1, cs('m2', 1))
        self.assertArrayEqual(cs1, cs('m2', 2))
        self.assertArrayEqual(cs1, cs('m2', 3))

    def test_valence(self):
        cmd.fragment('gly')
        cmd.valence(2, 'ID 0', 'ID 1')
        self.assertEqual(2, cmd.get_model('ID 0+1').bond[0].order)

    def test_vdw_fit(self):
        cmd.vdw_fit
        self.skipTest("TODO")

    @testing.requires_version('1.8.3.1')
    def test_set_discrete(self):
        import pymol

        cmd.fragment('ala', 'm1')
        cmd.create('m1', 'm1', 1, 2)

        self.assertEqual(0, cmd.count_discrete('*'))
        self.assertEqual(2, cmd.count_states())
        self.assertEqual(10, cmd.count_atoms())

        pymol.editing.set_discrete('m1', 1)

        self.assertEqual(1, cmd.count_discrete('*'))
        self.assertEqual(2, cmd.count_states())
        self.assertEqual(20, cmd.count_atoms())

        pymol.editing.set_discrete('m1', 0)

        self.assertEqual(0, cmd.count_discrete('*'))
        self.assertEqual(2, cmd.count_states())
        self.assertEqual(10, cmd.count_atoms())
