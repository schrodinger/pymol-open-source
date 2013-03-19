
'''
Testing atom properties for simple getting/setting and loading from sdf and mae files
'''

import unittest
from pymol import cmd, testing, stored

@testing.requires('properties')
class TestCopyProperties(testing.PyMOLTestCase):

    def testCopyObjectProperties(self):
        cmd.load(self.datafile('1molecule.mae'), 'test', load_properties='*')
        objs = cmd.get_object_list()
        for obj in objs:
            obj_copy = '%s_copy' % obj
            cmd.create(obj_copy, obj, copy_properties=True)
            props1 = cmd.get_property_list(obj)
            props2 = cmd.get_property_list(obj_copy)
            self.assertEqual(set(props1), set(props2))
            for prop in props1:
                prop1 = cmd.get_property(prop, obj)
                prop2 = cmd.get_property(prop, obj_copy)
                self.assertEqual(prop1, prop2)

    def testCopyAtomProperties(self):
        cmd.load(self.datafile('1molecule.mae'), 'test', load_properties='*', load_atom_properties='*')
        objs = cmd.get_object_list()
        for obj in objs:
            obj_copy = '%s_copy' % obj
            cmd.create(obj_copy, obj, copy_properties=True)

            stored.proplookup1 = {}
            cmd.iterate(obj,  'stored.proplookup1[index] = property')
            stored.proplookup2 = {}
            cmd.iterate(obj_copy,  'stored.proplookup2[index] = property')
            self.assertEqual(stored.proplookup1, stored.proplookup2)

