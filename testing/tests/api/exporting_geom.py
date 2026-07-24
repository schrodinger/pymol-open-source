'''
unit tests for pymol.exporting geometry formats
'''

import shutil
import struct
import subprocess
import zipfile

from pymol import cmd, testing

def file_get_contents(filename, mode='r'):
    with open(filename, mode) as handle:
        return handle.read()

class TestExportingGeom(testing.PyMOLTestCase):

    def testVRML(self):
        cmd.fragment('gly')
        for rep in ['spheres', 'sticks', 'surface']:
            cmd.show_as(rep)
            with testing.mktemp('.wrl') as filename:
                cmd.save(filename)
                contents = file_get_contents(filename)
                self.assertTrue(contents.startswith('#VRML V2'))

    @testing.requires_version('1.8')
    def testCOLLADA(self):
        cmd.fragment('gly')
        for rep in ['spheres', 'sticks', 'surface']:
            cmd.show_as(rep)
            with testing.mktemp('.dae') as filename:
                cmd.save(filename)
                contents = file_get_contents(filename)
                self.assertTrue('<COLLADA' in contents)

    def testUSDA(self):
        cmd.fragment('gly')
        for rep, schema in [
                ('spheres', 'def Sphere'),
                ('sticks', 'def Cylinder'),
                ('surface', 'def Mesh')]:
            cmd.show_as(rep)
            with testing.mktemp('.usda') as filename:
                cmd.save(filename)
                contents = file_get_contents(filename)
                self.assertTrue(contents.startswith('#usda 1.0'))
                self.assertIn('metersPerUnit = 1e-10', contents)
                self.assertIn(schema, contents)

    def testUSDZ(self):
        if shutil.which('usdcat') is None:
            self.skipTest('usdcat not available')

        cmd.fragment('gly')
        usdchecker = shutil.which('usdchecker')

        for rep in ['spheres', 'sticks', 'surface']:
            cmd.show_as(rep)

            with testing.mktemp('.usdz') as filename:
                cmd.save(filename)

                with zipfile.ZipFile(filename) as archive:
                    self.assertEqual(archive.namelist(), ['scene.usdc'])
                    info = archive.getinfo('scene.usdc')
                    self.assertEqual(info.compress_type, zipfile.ZIP_STORED)

                    with open(filename, 'rb') as handle:
                        handle.seek(info.header_offset)
                        header = handle.read(30)
                    fields = struct.unpack('<IHHHHHIIIHH', header)
                    data_offset = (
                        info.header_offset + 30 + fields[-2] + fields[-1])
                    self.assertEqual(data_offset % 64, 0)

                if usdchecker:
                    result = subprocess.run(
                        [usdchecker, '--arkit', filename],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True)
                    self.assertEqual(result.returncode, 0, result.stdout)

    @testing.requires('incentive')
    @testing.requires_version('2.1')
    def testSTL(self):
        cmd.fragment('gly')
        for rep in ['spheres', 'sticks', 'surface']:
            cmd.show_as(rep)
            with testing.mktemp('.stl') as filename:
                cmd.save(filename)
                contents = file_get_contents(filename, 'rb')
                # 80 bytes header
                # 4 bytes (uint32) number of triangles
                self.assertTrue(len(contents) > 84)
