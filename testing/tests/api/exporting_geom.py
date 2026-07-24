'''
unit tests for pymol.exporting geometry formats
'''

import re
import shutil
import struct
import subprocess
import unittest.mock
import zipfile

import pymol
from pymol import cmd, testing

def file_get_contents(filename, mode='r'):
    with open(filename, mode) as handle:
        return handle.read()

# two atoms with anisotropic displacement parameters, which are drawn as
# ellipsoids by the "ellipsoids" representation
v_pdbstr_anisou = (
    'ATOM      1  N   GLU A 114      24.832  -7.270  -5.728  1.00 33.91           N  \n'
    'ANISOU    1  N   GLU A 114     6968   3709   2207   -518    495    146       N  \n'
    'ATOM      2  CA  GLU A 114      25.839  -6.416  -5.102  1.00 33.68           C  \n'
    'ANISOU    2  CA  GLU A 114     6957   3558   2282   -477    534    166       C  \n'
    'END\n')

def usda_translations(contents):
    '''
    Centers of all prims which use a translate-only transform
    '''
    return sorted(
        tuple(round(float(value), 3) for value in match.split(','))
        for match in re.findall(
            r'float3 xformOp:translate = \(([^)]*)\)', contents))

def usda_matrix_translations(contents):
    '''
    Translation rows of all prims which use a full matrix transform
    '''
    pattern = re.compile(
        r'matrix4d xformOp:transform = \('
        r'\s*\([^)]*\),\s*\([^)]*\),\s*\([^)]*\),\s*\(([^)]*)\)')
    return sorted(
        tuple(round(float(value), 3) for value in match.split(',')[:3])
        for match in pattern.findall(contents))

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

    @testing.foreach(0, 1)
    def testUSDAEllipsoid(self, geometry_export_mode):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1', zoom=0)
        cmd.show_as('spheres')
        cmd.show('ellipsoids')
        cmd.turn('x', 35)
        cmd.turn('y', 50)
        cmd.set('geometry_export_mode', geometry_export_mode)

        with testing.mktemp('.usda') as filename:
            cmd.save(filename)
            contents = file_get_contents(filename)

        # every ellipsoid is centered on the atom which also carries a sphere,
        # so both must be written in the same coordinate space
        ellipsoids = usda_matrix_translations(contents)
        self.assertEqual(len(ellipsoids), 2)
        self.assertEqual(usda_translations(contents), ellipsoids)

        model = sorted(
            tuple(round(float(value), 3) for value in coord)
            for coord in cmd.get_coords('m1'))

        if geometry_export_mode:
            self.assertEqual(ellipsoids, model)
        else:
            self.assertNotEqual(ellipsoids, model)

    def testUSDZTimeout(self):
        from pymol import exporting

        expired = subprocess.TimeoutExpired('usdcat', 0.25)

        with unittest.mock.patch.object(
                exporting, '_find_usdcat', lambda: 'usdcat'), \
                unittest.mock.patch.object(
                    subprocess, 'run', side_effect=expired):
            with self.assertRaises(pymol.CmdException) as caught:
                exporting._convert_usda_to_usdc('#usda 1.0\n', timeout=0.25)

        self.assertIn('timed out', str(caught.exception))

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
