
'''
Testing settings for simple getting/setting 
'''

import unittest
from pymol import cmd, testing, stored

@testing.requires('properties')
class TestSettings3f(testing.PyMOLTestCase):
    @testing.foreach.product((tuple, list), (((1.,2.,3.), (3.,4.,5.),), ((6.,7.,8.),(9.,0.,1.),),),)
    def test_set_object_settings_3f(self, f, data):
        m1 = "pseudo01"
        m2 = "pseudo02"
        cmd.pseudoatom(m1)
        cmd.pseudoatom(m2)
        lp1, lp2 = lp = map(f, data)
        cmd.set("label_position", lp1, m1)
        cmd.set("label_position", lp2, m2)
        stored.pos = {}
        cmd.iterate("all", "stored.pos[model]=s['label_position']")
        self.assertEqual(tuple(lp1), stored.pos[m1])
        self.assertEqual(tuple(lp2), stored.pos[m2])

    @testing.foreach.product((tuple, list), (((1.,2.,3.), (3.,4.,5.),), ((6.,7.,8.),(9.,0.,1.),),),)
    def test_set_both_object_and_global_settings_3f(self, f, data):
        m1 = "pseudo01"
        m2 = "pseudo02"
        m3 = "pseudo03"
        lp1, lp2 = lp = map(f, data)
        cmd.pseudoatom(m1)
        cmd.pseudoatom(m2)
        cmd.pseudoatom(m3)
        cmd.set("label_position", lp1)
        cmd.set("label_position", lp2, m2)
        stored.pos = {}
        cmd.iterate("all", "stored.pos[model]=s['label_position']")
        self.assertEqual(tuple(lp1), stored.pos[m1])
        self.assertEqual(tuple(lp2), stored.pos[m2])
        self.assertEqual(tuple(lp1), stored.pos[m3])

    @testing.foreach.product((tuple, list), (((1.,2.,3.), (3.,4.,5.), (6.,7.,8.),), ((6.,7.,8.),(9.,0.,1.), (2.,3.,4.),),),)
    def test_set_atom_object_and_global_settings_3f(self, f, data):
        m1 = "pseudo01"
        m2 = "pseudo02"
        m3 = "pseudo03"
        lp1, lp2, lp3 = lp = map(f, data)
        cmd.pseudoatom(m1)
        cmd.pseudoatom(m1)
        cmd.pseudoatom(m2)
        cmd.pseudoatom(m3)
        cmd.set("label_position", lp1, "(index 2)")
        cmd.set("label_position", lp2)
        cmd.set("label_position", lp3, m2)
        stored.pos = {}
        cmd.iterate("all", "stored.pos['%s-%s' % (model, index)]=s['label_position']")
        self.assertEqual(tuple(lp1), stored.pos['%s-%s' % (m1,2)])
        self.assertEqual(tuple(lp2), stored.pos['%s-%s' % (m1,1)])
        self.assertEqual(tuple(lp3), stored.pos['%s-%s' % (m2,1)])
        self.assertEqual(tuple(lp2), stored.pos['%s-%s' % (m3,1)])

    @testing.foreach.product((tuple, list), (((1.,2.,3.), (3.,4.,5.),), ((6.,7.,8.),(9.,0.,1.),),),)
    def test_atom_state_settings_3f(self, f, data):
        m1 = "pseudo01"
        lp1, lp2 = lp = map(f, data)
        cmd.pseudoatom(m1)
        cmd.pseudoatom(m1)
        cmd.create(m1, m1, 1, 2)
        stored.pos = {}
        stored.lp1 = lp1
        stored.lp2 = lp2
        stored.origp = None
        # get origp (should be 0.,0.,0.
        cmd.alter("(index 1)", "stored.origp = s['label_placement_offset']")
        # change atom-state setting to lp1 for atom in both states
        cmd.alter_state(0, "(index 1)", "s['label_placement_offset']=stored.lp1")
        # change atom-level setting for all atoms to something else
        cmd.alter("all", "s['label_placement_offset']=stored.lp2")
        # get atom-state settings
        cmd.iterate_state(0, "all", "stored.pos['%s-%s-%s' % (model, state, index)]=s['label_placement_offset']")
        # atom-state setting should be lp1 for atom 1
        self.assertEqual(tuple(lp1), stored.pos['%s-%s-%s' % (m1, 1, 1)])
        self.assertEqual(tuple(lp1), stored.pos['%s-%s-%s' % (m1, 2, 1)])
        # atom setting should override to lp2 for atom 2 in both states
        self.assertEqual(tuple(lp2), stored.pos['%s-%s-%s' % (m1, 1, 2)])
        self.assertEqual(tuple(lp2), stored.pos['%s-%s-%s' % (m1, 2, 2)])
        # unset all atom-level settings, atom-state settings should still exist
        cmd.alter("all", "s['label_placement_offset']=None")
        stored.pos = {}
        cmd.iterate_state(0, "all", "stored.pos['%s-%s-%s' % (model, state, index)]=s['label_placement_offset']")
        self.assertEqual(tuple(lp1), stored.pos['%s-%s-%s' % (m1, 1, 1)])
        self.assertEqual(tuple(lp1), stored.pos['%s-%s-%s' % (m1, 2, 1)])
        self.assertEqual(tuple(stored.origp), stored.pos['%s-%s-%s' % (m1, 1, 2)])
        self.assertEqual(tuple(stored.origp), stored.pos['%s-%s-%s' % (m1, 2, 2)])

    @testing.foreach.product(((list),), ((1.,2.,3.), (3.,4.,5.),))
    def test_unsetting_atom_setting(self, f, data):
        m1 = "pseudo01"
        lp = f(data)
        cmd.pseudoatom(m1)
        stored.lp = lp
        stored.pos = {}
        cmd.iterate("all", "stored.pos['%s-%s' % (model, index)] = s['label_placement_offset']")

        cmd.alter("all", "s['label_placement_offset']=stored.lp")
        stored.pos2 = {}
        cmd.alter("all", "s['label_placement_offset']=None")
        cmd.iterate("all", "stored.pos2['%s-%s' % (model, index)] = s['label_placement_offset']")
        self.assertEqual(stored.pos, stored.pos2)
