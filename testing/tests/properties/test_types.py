
'''
Testing setting/getting properties for different types
'''

import random
import unittest
from pymol import cmd, testing

try:
    from pymol.properties import PROPERTY_AUTO, PROPERTY_BOOLEAN, \
     PROPERTY_INT, PROPERTY_FLOAT, PROPERTY_STRING, PROPERTY_COLOR
except ImportError:
    pass

@testing.requires('properties')
class TestProperties(testing.PyMOLTestCase):
    def testBooleanType(self):
        cmd.fragment('ala')
        cmd.set_property('bool_property_1', True, 'ala')
        cmd.set_property('bool_property_2', False, 'ala')
        self.assertEqual(cmd.get_property('bool_property_1', 'ala'), True)
        self.assertEqual(cmd.get_property('bool_property_2', 'ala'), False)

        # test setting boolean type by string
        cmd.set_property('bool_property_1', "True", 'ala', proptype=PROPERTY_BOOLEAN)
        cmd.set_property('bool_property_2', "True", 'ala', proptype=PROPERTY_BOOLEAN)
        cmd.set_property('bool_property_3', "True", 'ala', proptype=PROPERTY_BOOLEAN)
        cmd.set_property('bool_property_4', "False", 'ala', proptype=PROPERTY_BOOLEAN)
        self.assertEqual(cmd.get_property('bool_property_1', 'ala'), True)
        self.assertEqual(cmd.get_property('bool_property_2', 'ala'), True)   # overwritting 
        self.assertEqual(cmd.get_property('bool_property_3', 'ala'), True)
        self.assertEqual(cmd.get_property('bool_property_4', 'ala'), False)

    @testing.foreach((random.randint(1, 100000), random.randint(1, 100000)), (-random.randint(1, 100000), -random.randint(1, 100000)))
    def testIntegerType(self, val1, val2):
        cmd.fragment('ala')
        cmd.set_property('int_property_1', val1, 'ala')
        cmd.set_property('int_property_2', val2, 'ala')
        self.assertEqual(cmd.get_property('int_property_1', 'ala'), val1)
        self.assertEqual(cmd.get_property('int_property_2', 'ala'), val2)

        # test setting boolean type by string
        cmd.set_property('int_property_1', str(val1), 'ala', proptype=PROPERTY_INT)
        cmd.set_property('int_property_2', str(val1), 'ala', proptype=PROPERTY_INT)
        cmd.set_property('int_property_3', str(val1), 'ala', proptype=PROPERTY_INT)
        cmd.set_property('int_property_4', str(val2), 'ala', proptype=PROPERTY_INT)
        self.assertEqual(cmd.get_property('int_property_1', 'ala'), val1)
        self.assertEqual(cmd.get_property('int_property_2', 'ala'), val1)   # overwritting 
        self.assertEqual(cmd.get_property('int_property_3', 'ala'), val1)
        self.assertEqual(cmd.get_property('int_property_4', 'ala'), val2)

    @testing.foreach((random.random() * 100000, random.random() * 100000), (-random.random() * 100000, -random.random() * 100000))
    def testFloatType(self, val1, val2):
        cmd.fragment('ala')
        cmd.set_property('float_property_1', val1, 'ala')
        cmd.set_property('float_property_2', val2, 'ala')
        self.assertAlmostEqual(cmd.get_property('float_property_1', 'ala'), val1)
        self.assertAlmostEqual(cmd.get_property('float_property_2', 'ala'), val2)

        # test setting boolean type by string
        cmd.set_property('float_property_1', str(val1), 'ala', proptype=PROPERTY_FLOAT)
        cmd.set_property('float_property_2', str(val1), 'ala', proptype=PROPERTY_FLOAT)
        cmd.set_property('float_property_3', str(val1), 'ala', proptype=PROPERTY_FLOAT)
        cmd.set_property('float_property_4', str(val2), 'ala', proptype=PROPERTY_FLOAT)
        self.assertAlmostEqual(cmd.get_property('float_property_1', 'ala'), val1)
        self.assertAlmostEqual(cmd.get_property('float_property_2', 'ala'), val1)   # overwritting 
        self.assertAlmostEqual(cmd.get_property('float_property_3', 'ala'), val1)
        self.assertAlmostEqual(cmd.get_property('float_property_4', 'ala'), val2)

    @testing.foreach(('test string 1', 'test string 2'))
    def testStringType(self, str1, str2):
        cmd.fragment('ala')
        cmd.set_property('string_property_1', str1, 'ala')
        cmd.set_property('string_property_2', str2, 'ala')
        self.assertEqual(cmd.get_property('string_property_1', 'ala'), str1)
        self.assertEqual(cmd.get_property('string_property_2', 'ala'), str2)

        # test setting boolean type by string
        cmd.set_property('string_property_1', str1, 'ala', proptype=PROPERTY_STRING)
        cmd.set_property('string_property_2', str1, 'ala', proptype=PROPERTY_STRING)
        self.assertEqual(cmd.get_property('string_property_1', 'ala'), str1)
        self.assertEqual(cmd.get_property('string_property_2', 'ala'), str1)   # overwritting 

    @testing.foreach(('red', 'blue'), ('0xff0000', '0x0000ff'), ([ 1., 0., 0.], [0., 0., 1.]), (4,5))
    def testColorType(self, color1, color2):
        # colors are saved as integers, but can be accessed/set by either strings or hex
        from pymol.querying import get_color_index_from_string_or_list
        val1 = get_color_index_from_string_or_list(color1, _self=cmd)
        val2 = get_color_index_from_string_or_list(color2, _self=cmd)
        cmd.fragment('ala')
        cmd.set_property('color_property_1', color1, 'ala', proptype=PROPERTY_COLOR)
        cmd.set_property('color_property_2', color2, 'ala', proptype=PROPERTY_COLOR)
        self.assertColorEqual(cmd.get_property('color_property_1', 'ala'), val1)
        self.assertColorEqual(cmd.get_property('color_property_2', 'ala'), val2)
        self.assertColorEqual(cmd.get_property('color_property_1', 'ala'), color1)
        self.assertColorEqual(cmd.get_property('color_property_2', 'ala'), color2)
        cmd.set_property('color_property_1', str(val1), 'ala', proptype=PROPERTY_COLOR)
        cmd.set_property('color_property_2', str(val1), 'ala', proptype=PROPERTY_COLOR)
        cmd.set_property('color_property_3', str(val1), 'ala', proptype=PROPERTY_COLOR)
        cmd.set_property('color_property_4', str(val2), 'ala', proptype=PROPERTY_COLOR)
        self.assertColorEqual(cmd.get_property('color_property_1', 'ala'), color1)
        self.assertColorEqual(cmd.get_property('color_property_1', 'ala'), val1)
        self.assertColorEqual(cmd.get_property('color_property_2', 'ala'), color1)   # overwritting 
        self.assertColorEqual(cmd.get_property('color_property_2', 'ala'), val1)
        self.assertColorEqual(cmd.get_property('color_property_3', 'ala'), color1)
        self.assertColorEqual(cmd.get_property('color_property_3', 'ala'), val1)
        self.assertColorEqual(cmd.get_property('color_property_4', 'ala'), color2)
        self.assertColorEqual(cmd.get_property('color_property_4', 'ala'), val2)

    def testColorTypeCmd(self):
        cmd.fragment('ala')
        cmd.do('set_property color_prop, red, ala, proptype=5', echo=0)
        prop = cmd.get_property('color_prop', 'ala')
        self.assertColorEqual(prop, 'red')

        
        
