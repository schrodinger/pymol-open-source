'''
Testing: pymol.creating
'''

import pymol
from pymol import cmd, testing, stored

import unittest

class TestCreating(testing.PyMOLTestCase):

    def testGroup(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3')
        m_list = []
        cmd.group('g1', 'm1 m2')
        cmd.iterate('g1', 'm_list.append(model)', space=locals())
        self.assertItemsEqual(m_list, ['m1', 'm2'])
        m_list = []
        cmd.ungroup('m2')
        cmd.iterate('g1', 'm_list.append(model)', space=locals())
        self.assertItemsEqual(m_list, ['m1'])

    def testUngroup(self):
        # see testGroup
        pass

    @testing.requires_version('1.7.7')
    def testGroupRename(self):
        '''
        Rename group entries which have a group name prefix when renaming the group
        '''
        cmd.set('group_auto_mode', 2)
        cmd.pseudoatom('g1.m1')
        cmd.pseudoatom('g1.m2')
        cmd.set_name('g1', 'g2')
        m_list = []
        cmd.iterate('g2', 'm_list.append(model)', space=locals())
        self.assertItemsEqual(m_list, ['g2.m1', 'g2.m2'])

    def testIsomesh(self):
        cmd.viewport(100, 100)

        cmd.fragment('gly', 'm1')
        cmd.set('gaussian_b_floor', 30)
        cmd.set('mesh_width', 5)
        cmd.map_new('map')
        cmd.delete('m1')

        # make mesh
        cmd.isomesh('mesh', 'map')
        cmd.isomesh('mesh', 'map', source_state=1, state=-2)

        ## check mesh presence by color
        meshcolor = 'red'
        cmd.color(meshcolor, 'mesh')
        self.ambientOnly()
        self.assertImageHasColor(meshcolor)
        cmd.frame(2)
        self.assertEqual(cmd.get_state(), 2)
        self.assertImageHasColor(meshcolor)

        with testing.mktemp('.pse') as filename:
            cmd.save(filename)

            cmd.delete('*')
            self.assertImageHasNotColor(meshcolor)

            cmd.load(filename)
            self.assertImageHasColor(meshcolor)
            cmd.frame(2)
            self.assertEqual(cmd.get_state(), 2)
            self.assertImageHasColor(meshcolor)

    def testIsosurface(self):
        cmd.viewport(100, 100)

        cmd.fragment('gly', 'm1')
        cmd.set('gaussian_b_floor', 30)
        cmd.set('mesh_width', 5)
        cmd.map_new('map')
        cmd.delete('m1')

        # make mesh
        cmd.isosurface('surface', 'map')
        cmd.isosurface('surface', 'map', source_state=1, state=-2)

        ## check mesh presence by color
        meshcolor = 'red'
        cmd.color(meshcolor, 'surface')
        self.ambientOnly()
        self.assertImageHasColor(meshcolor)
        cmd.frame(2)
        self.assertEqual(cmd.get_state(), 2)
        self.assertImageHasColor(meshcolor)

        with testing.mktemp('.pse') as filename:
            cmd.save(filename)

            cmd.delete('*')
            self.assertImageHasNotColor(meshcolor)

            cmd.load(filename)
            self.assertImageHasColor(meshcolor)
            cmd.frame(2)
            self.assertEqual(cmd.get_state(), 2)
            self.assertImageHasColor(meshcolor)

    def testIsodot(self):
        cmd.viewport(100, 100)

        cmd.fragment('gly', 'm1')
        cmd.set('gaussian_b_floor', 30)
        cmd.set('mesh_width', 5)
        cmd.map_new('map')
        cmd.delete('m1')

        # make mesh
        cmd.isodot('dot', 'map')
        cmd.isodot('dot', 'map', source_state=1, state=-2)

        ## check mesh presence by color
        meshcolor = 'red'
        cmd.color(meshcolor, 'dot')
        self.ambientOnly()
        self.assertImageHasColor(meshcolor)
        cmd.frame(2)
        self.assertEqual(cmd.get_state(), 2)
        self.assertImageHasColor(meshcolor)

        with testing.mktemp('.pse') as filename:
            cmd.save(filename)

            cmd.delete('*')
            self.assertImageHasNotColor(meshcolor)

            cmd.load(filename)
            self.assertImageHasColor(meshcolor)
            cmd.frame(2)
            self.assertEqual(cmd.get_state(), 2)
            self.assertImageHasColor(meshcolor)

    def testIsolevel(self):
        cmd.viewport(100, 100)

        cmd.fragment('gly', 'm1')
        cmd.set('gaussian_b_floor', 30)
        cmd.set('mesh_width', 5)
        cmd.map_new('map')
        cmd.delete('m1')

        # make mesh
        cmd.isodot('dot', 'map')
        cmd.isodot('dot', 'map', source_state=1, state=-2)

        ## check mesh presence by color
        meshcolor = 'red'
        cmd.color(meshcolor, 'dot')
        self.ambientOnly()
        self.assertImageHasColor(meshcolor)
        cmd.frame(2)
        self.assertEqual(cmd.get_state(), 2)
        self.assertImageHasColor(meshcolor)

        for contourlvl in range(7):
            cmd.isolevel('dot', contourlvl)
            self.assertImageHasColor(meshcolor)

        cmd.isolevel('dot', 10)
        self.assertImageHasNotColor(meshcolor)

    @unittest.skip("Not yet implemented.")
    def testGradient(self):
        pass

    def testCopy(self):
        cmd.fragment('ala', 'm1')
        cmd.copy('m2', 'm1')
        self.assertEqual(
                cmd.count_atoms('m1'),
                cmd.count_atoms('m2'))

    @testing.requires_version('2.4')
    def testCopyMap(self):
        cmd.load(self.datafile('emd_1155.ccp4'), 'map1')
        cmd.copy('map2', 'map1')
        self.assertEqual(cmd.get_symmetry('map1'), cmd.get_symmetry('map2'))
        self.assertArrayEqual(cmd.get_volume_field('map1'), cmd.get_volume_field('map2'))

    def testSymexp(self):
        cmd.load(self.datafile('1oky.pdb.gz'), 'm1')
        n = cmd.count_atoms()
        cmd.symexp('s', 'm1', '%m1 & resi 283', 20.0)
        x = cmd.get_object_list()
        self.assertEqual(x, [
            'm1',
            's01000000',
            's03000000',
            's04000000',
            ])
        self.assertEqual(n * 4, cmd.count_atoms())

    @unittest.skip("Not yet implemented.")
    def testFragment(self):
        pass

    def testCreate(self):
        cmd.fragment("ala")
        cmd.create("foo", "ala", 1, 1)

        names = cmd.get_names()
        names.sort()

        self.assertEquals(names, ["ala", "foo"])
        
        # TODO: Create a function to compare to chempy models
        # m1 = cmd.get_model("foo")
        # m2 = cmd.get_model("ala")
        # self.assertChempyModelEquals(m1, m2)

        self.assertEquals(cmd.fit("ala", "foo"), 0.0)

        self.assertEquals(cmd.count_states(), 1)

    def testCreateMultipleStates(self):
        cmd.fragment("ala")

        # simple creation

        self.assertEquals(cmd.count_states(), 1)

        cmd.create("foo", "ala", 1, 1)

        self.assertEquals(cmd.count_states(), 1)

        # creation into state 1

        cmd.create("foo", "ala", 1, 1)

        self.assertEquals(cmd.count_states(), 1)

        # creation into state 2

        cmd.create("foo", "ala", 1, 2)

        self.assertEquals(cmd.count_states(), 2)
        
    def create_many_states(self, fragment, obj, count):
        cmd.fragment(fragment, obj)
        for x in range(2, count + 1):
            cmd.create(obj, obj, 1, x)

    def testCreateMany(self):

        # 10 states
        self.create_many_states("thr", "foo", 10)

        self.assertEquals(cmd.count_states("foo"), 10)
        self.assertEquals(cmd.count_states(), 10)

        # 100 states
        self.create_many_states("trp", "bar", 100)

        self.assertEquals(cmd.count_states("foo"), 10)
        self.assertEquals(cmd.count_states("bar"), 100)
        self.assertEquals(cmd.count_states(), 100)

        # 1000 states
        
        self.create_many_states("ile", "bam", 1000)

        self.assertEquals(cmd.count_states(), 1000)
        
    @testing.requires_version('1.6')
    def testCreateTargetState(self):
        self.create_many_states("gly", "m1", 4)
        cmd.frame(3)

        cmd.create('m2', 'm1', 0, 0)
        self.assertEqual(cmd.count_states('m2'), 4)
        cmd.delete('m2')

        cmd.create('m2', 'm1', 1, 0)
        self.assertEqual(cmd.count_states('m2'), 1)
        cmd.delete('m2')

        cmd.create('m2', 'm1', 2, 0)
        self.assertEqual(cmd.count_states('m2'), 2)
        cmd.delete('m2')

        cmd.create('m2', 'm1', 0, 2)
        self.assertEqual(cmd.count_states('m2'), 5)
        cmd.delete('m2')

        cmd.create('m2', 'm1', 1, 2)
        self.assertEqual(cmd.count_states('m2'), 2)
        cmd.delete('m2')

        cmd.create('m2', 'm1', 2, 2)
        self.assertEqual(cmd.count_states('m2'), 2)
        cmd.delete('m2')

        cmd.create('m2', 'm1', 0, -1)
        self.assertEqual(cmd.count_states('m2'), 4)
        cmd.create('m2', 'm1', 0, -1)
        self.assertEqual(cmd.count_states('m2'), 8)
        cmd.create('m2', 'm1', 2, -1)
        self.assertEqual(cmd.count_states('m2'), 9)

    @testing.requires_version('2.1')
    def testCreateCurrentSourceState(self):
        self.create_many_states("gly", "m1", 4)
        cmd.frame(3)

        cmd.create('m2', 'm1', -1, -1)
        self.assertEqual(cmd.count_states('m2'), 1)
        cmd.delete('m2')

        cmd.create('m2', 'm1', -1, 0)
        self.assertEqual(cmd.count_states('m2'), 3)
        cmd.delete('m2')

        cmd.create('m2', 'm1', -1, 1)
        self.assertEqual(cmd.count_states('m2'), 1)
        cmd.delete('m2')

        cmd.create('m2', 'm1', -1, 2)
        self.assertEqual(cmd.count_states('m2'), 2)
        cmd.delete('m2')

    # PYMOL-3372
    @testing.foreach(0, 1)
    @testing.requires_version('2.4')
    def testCreateDiscrete(self, discrete_source):
        self.create_many_states("methanol", "m1", 2)
        self.create_many_states("benzene", "m2", 2)
        pymol.editing.set_discrete("m1", discrete=discrete_source)
        pymol.editing.set_discrete("m2", discrete=discrete_source)

        # create new
        cmd.create('dm', 'm1', 0, -1, discrete=1)
        self.assertTrue(cmd.count_discrete('dm'))
        self.assertEqual(cmd.count_states('dm'), 2)
        self.assertEqual(cmd.count_atoms('dm'), 2 * 6)
        self.assertEqual(len(cmd.get_bonds('dm')), 5)

        # append to existing
        cmd.create('dm', 'm2', 0, -1, discrete=1)
        self.assertTrue(cmd.count_discrete('dm'))
        self.assertEqual(cmd.count_states('dm'), 4)
        self.assertEqual(cmd.count_atoms('dm'), 2 * 6 + 2 * 12)
        self.assertEqual(len(cmd.get_bonds('dm', 1)), 5)
        self.assertEqual(len(cmd.get_bonds('dm', 4)), 12)
        self.assertEqual(len(cmd.get_bonds('dm', 0)), 2 * 5 + 2 * 12)

    def testExtract(self):
        cmd.fragment('ala', 'm1')
        cmd.extract('m2', 'elem C')
        self.assertEqual(cmd.count_atoms('m1'), 7)
        self.assertEqual(cmd.count_atoms('m2'), 3)
    
    @unittest.skip("Not yet implemented.")
    def testUnquote(self):
        pass

    def testPseudoatom(self):
        cmd.set('retain_order')

        cmd.pseudoatom('m1', '', 'A1', 'B2', 1, 'D',
                'E5', 'F', 7.0, 0, 9.0, 0.1, '11', 'foo12',
                (13, 14, 15), 1)
        cmd.pseudoatom('m2', '', 'A1', 'B2', 2, 'A',
                'E5', 'F', 7.0, 0, 9.0, 0.1, '12', 'bar12',
                (13, 10, 15), 1)

        # append to m1, pos and vdw vom selection
        cmd.pseudoatom('m1', 'all', 'A2', 'B2', 3, 'D',
                'E5', 'F', -1, 1, 19.0, 0.2, '13', 'com12',
                None, 1, 'extent')

        r_list = []
        f_indices = [7, 9, 10, 13, 14, 15]

        cmd.iterate_state(1, 'all',
                'r_list.append([model, name, resn, resi, chain, segi, elem, '
                'vdw, type, b, q, color, label, x, y, z])', space=locals())

        for ref, values in zip(r_list, [
            ['m1', 'A1', 'B2', '1', 'D', 'E5', 'F', 7.0, 'ATOM', 9.0, 0.1, 11, 'foo12', 13.0, 14.0, 15.0],
            ['m1', 'A2', 'B2', '3', 'D', 'E5', 'F', 2.0, 'HETATM', 19.0, 0.2, 13, 'com12', 13.0, 12.0, 15.0],
            ['m2', 'A1', 'B2', '2', 'A', 'E5', 'F', 7.0, 'ATOM', 9.0, 0.1, 12, 'bar12', 13.0, 10.0, 15.0],
            ]):
            for i in f_indices:
                values[i] = round(values[i], 1)
                ref[i] = round(ref[i], 1)
            self.assertEqual(ref, values)

    def testPseudoatomName(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m1')
        self.assertEqual(2, cmd.count_atoms())
        self.assertEqual(1, cmd.count_atoms('name PS1'))
        self.assertEqual(1, cmd.count_atoms('name PS2'))

    @testing.requires_version('2.3')
    def test_set_raw_alignment(self):
        cmd.fab('ACDEF', 'm1')
        cmd.fab('CDE', 'm2')
        index_m1 = [('m1', 12), ('m1', 23), ('m1', 35)]
        index_m2 = [('m2',  2), ('m2', 13), ('m2', 25)]
        raw = [list(t) for t in zip(index_m1, index_m2)]
        cmd.set_raw_alignment('aln', raw)
        self.assertEqual(cmd.index('m1 & aln'), index_m1)
        self.assertEqual(cmd.index('m2 & aln'), index_m2)
