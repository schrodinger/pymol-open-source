'''
unit tests for pymol.exporting geometry formats
'''

import math
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

# an ANISOU tensor which is not positive definite (eigenvalues -3000, 6000 and
# 15000), so that one ellipsoid semi-axis collapses to zero
v_pdbstr_anisou_indefinite = (
    'ATOM      1  N   GLU A 114      24.832  -7.270  -5.728  1.00 33.91           N  \n'
    'ANISOU    1  N   GLU A 114     6000   6000   6000   9000      0      0       N  \n'
    'END\n')

def usda_translations(contents):
    '''
    Centers of all prims which use a translate-only transform
    '''
    return sorted(
        tuple(round(float(value), 3) for value in match.split(','))
        for match in re.findall(
            r'float3 xformOp:translate = \(([^)]*)\)', contents))

def usda_matrices(contents):
    '''
    Rows of all "matrix4d xformOp:transform" values
    '''
    pattern = re.compile(
        r'matrix4d xformOp:transform = \('
        r'\s*\(([^)]*)\),\s*\(([^)]*)\),\s*\(([^)]*)\),\s*\(([^)]*)\)')
    return [
        [tuple(float(value) for value in row.split(',')) for row in match]
        for match in pattern.findall(contents)]

def usda_matrix_translations(contents):
    '''
    Translation rows of all prims which use a full matrix transform
    '''
    return sorted(
        tuple(round(value, 3) for value in matrix[3][:3])
        for matrix in usda_matrices(contents))

def matrix_determinant(matrix):
    (a, b, c), (d, e, f), (g, h, i) = (row[:3] for row in matrix[:3])
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

def matrix_axis_lengths(matrix):
    return [math.sqrt(sum(value * value for value in row[:3]))
            for row in matrix[:3]]

def usda_prim_block(contents, name):
    '''
    Body of the prim with the given name, which must not have children
    '''
    match = re.search(
        r'def \w+ "' + re.escape(name) + r'"[^{]*\{(.*?)\n    \}', contents,
        re.DOTALL)
    return match and match.group(1)

def usda_vec3_array(block, declaration):
    match = re.search(
        re.escape(declaration) + r'\s*=\s*\[(.*?)\]', block, re.DOTALL)
    return [tuple(float(value) for value in item.split(','))
            for item in re.findall(r'\(([^)]*)\)', match.group(1))]

def usda_int_array(block, declaration):
    match = re.search(re.escape(declaration) + r'\s*=\s*\[([^\]]*)\]', block)
    return [int(value) for value in match.group(1).split(',')]

def ring_center_and_radius(points):
    center = [sum(point[axis] for point in points) / len(points)
              for axis in range(3)]
    radii = [math.dist(point, center) for point in points]
    return center, radii

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

        for matrix in usda_matrices(contents):
            self.assertRightHanded(matrix)

    def testUSDAEllipsoidMirrored(self):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1', zoom=0)
        cmd.show_as('ellipsoids')

        # a reflection turns the ellipsoid axes into a left-handed basis
        cmd.transform_object('m1', [
            -1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0])

        with testing.mktemp('.usda') as filename:
            cmd.save(filename)
            contents = file_get_contents(filename)

        matrices = usda_matrices(contents)
        self.assertEqual(len(matrices), 2)

        for matrix in matrices:
            self.assertRightHanded(matrix)

    def assertRightHanded(self, matrix):
        # renderers derive inward-pointing normals from a negative determinant
        lengths = matrix_axis_lengths(matrix)
        self.assertAlmostEqual(matrix_determinant(matrix),
                               lengths[0] * lengths[1] * lengths[2],
                               places=4)

    def testUSDAEllipsoidScale(self):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1', zoom=0)
        cmd.show_as('ellipsoids')

        lengths = {}

        for scale in (1.0, 2.0):
            cmd.set('ellipsoid_scale', scale)

            with testing.mktemp('.usda') as filename:
                cmd.save(filename)
                contents = file_get_contents(filename)

            matrices = usda_matrices(contents)
            self.assertEqual(len(matrices), 2)
            lengths[scale] = [
                sorted(matrix_axis_lengths(matrix)) for matrix in matrices]

        # semi-axes carry the ellipsoid size, they are not unit length
        for single, double in zip(lengths[1.0], lengths[2.0]):
            self.assertGreater(single[0], 0.0)
            for one, two in zip(single, double):
                self.assertAlmostEqual(two, one * 2.0, places=4)

    def testUSDAEllipsoidDegenerate(self):
        cmd.read_pdbstr(v_pdbstr_anisou_indefinite, 'm1', zoom=0)
        cmd.show_as('ellipsoids')

        with testing.mktemp('.usda') as filename:
            cmd.save(filename)
            contents = file_get_contents(filename)

        self.assertTrue(contents.startswith('#usda 1.0'))

        matrices = usda_matrices(contents)
        self.assertEqual(len(matrices), 1)

        for matrix in matrices:
            # a collapsed semi-axis must not make the transform singular
            lengths = matrix_axis_lengths(matrix)
            for length in lengths:
                self.assertTrue(math.isfinite(length))
                self.assertGreater(length, 0.0)
            self.assertGreater(matrix_determinant(matrix), 0.0)
            self.assertRightHanded(matrix)

            # the collapsed axis stays visually negligible
            self.assertAlmostEqual(min(lengths) / max(lengths), 1e-3, places=6)

        if shutil.which('usdcat') is None:
            self.skipTest('usdcat not available')

        # the layer must still be readable by OpenUSD
        with testing.mktemp('.usdz') as filename:
            cmd.save(filename)

    def testUSDAEllipsoidZeroScale(self):
        cmd.read_pdbstr(v_pdbstr_anisou, 'm1', zoom=0)
        cmd.show_as('ellipsoids')
        cmd.set('ellipsoid_scale', 0)

        with testing.mktemp('.usda') as filename:
            cmd.save(filename)
            contents = file_get_contents(filename)

        # ellipsoids without any extent are omitted rather than exported with
        # an all-zero (singular) transform
        self.assertTrue(contents.startswith('#usda 1.0'))
        self.assertEqual(usda_matrices(contents), [])
        self.assertNotIn('matrix4d', contents)

    def testUSDAAnalyticSolids(self):
        from pymol import cgo

        cmd.load_cgo([
            cgo.CONE, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0, 0.0,
            1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
            cgo.CYLINDER, 4.0, 0.0, 0.0, 4.0, 0.0, 5.0, 1.0,
            0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
        ], 'solids')

        with testing.mktemp('.usda') as filename:
            cmd.save(filename)
            contents = file_get_contents(filename)

        # cones and cylinders share one writer, so schema and prim name must
        # stay in sync
        solids = re.findall(r'def (Cone|Cylinder) "(\w+)_\d+"', contents)
        self.assertEqual(sorted(schema for schema, _ in solids),
                         ['Cone', 'Cylinder'])

        for schema, name in solids:
            self.assertEqual(schema, name)

    def save_cone_usda(self, r1, r2, c2, cap1=1.0, cap2=1.0):
        from pymol import cgo

        cmd.set('geometry_export_mode', 1)
        cmd.load_cgo([
            cgo.CONE, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, r1, r2,
            1.0, 0.0, 0.0, *c2, cap1, cap2,
        ], 'cone')

        with testing.mktemp('.usda') as filename:
            cmd.save(filename)
            return file_get_contents(filename)

    def testUSDAConeTwoColor(self):
        # a pointed cone with two endpoint colors cannot use UsdGeomCone,
        # which carries a single color
        contents = self.save_cone_usda(1.0, 0.0, (0.0, 0.0, 1.0))

        self.assertNotIn('def Cone', contents)

        block = usda_prim_block(contents, 'Cone_0')
        self.assertIsNotNone(block)

        colors = usda_vec3_array(block, 'color3f[] primvars:displayColor')
        points = usda_vec3_array(block, 'point3f[] points')
        self.assertEqual(len(colors), len(points))

        # the wide end keeps c1, the apex keeps c2
        by_z = {}
        for point, color in zip(points, colors):
            by_z.setdefault(round(point[2], 3), set()).add(color)

        levels = sorted(by_z)
        self.assertEqual(len(levels), 2)
        self.assertAlmostEqual(levels[1] - levels[0], 5.0, places=4)
        self.assertEqual(by_z[levels[0]], {(1.0, 0.0, 0.0)})
        self.assertEqual(by_z[levels[1]], {(0.0, 0.0, 1.0)})

    def testUSDAConeFrustum(self):
        # a truncated cone must keep both radii, not become one cylinder
        contents = self.save_cone_usda(2.0, 1.0, (1.0, 0.0, 0.0))

        self.assertNotIn('def Cylinder', contents)
        self.assertNotIn('def Cone', contents)

        block = usda_prim_block(contents, 'Cone_0')
        self.assertIsNotNone(block)

        points = usda_vec3_array(block, 'point3f[] points')
        counts = usda_int_array(block, 'int[] faceVertexCounts')
        indices = usda_int_array(block, 'int[] faceVertexIndices')

        self.assertEqual(len(indices), sum(counts))
        self.assertLess(max(indices), len(points))

        segments = counts.count(4)
        self.assertGreater(segments, 8)

        # both flat caps are triangle fans over the same segment count
        self.assertEqual(counts.count(3), 2 * segments)

        centers = []
        for offset, radius in [(0, 2.0), (segments, 1.0)]:
            center, radii = ring_center_and_radius(
                points[offset:offset + segments])
            centers.append(center)
            for value in radii:
                self.assertAlmostEqual(value, radius, places=4)

        self.assertAlmostEqual(math.dist(*centers), 5.0, places=4)

    def testUSDAConeUncapped(self):
        contents = self.save_cone_usda(
            2.0, 1.0, (1.0, 0.0, 0.0), cap1=0.0, cap2=0.0)

        block = usda_prim_block(contents, 'Cone_0')
        counts = usda_int_array(block, 'int[] faceVertexCounts')

        # only the lateral surface, no cap fans
        self.assertEqual(counts.count(3), 0)
        self.assertEqual(len(counts), counts.count(4))

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

    def convert_with_fake_usdcat(self, output):
        '''
        Run the usdcat conversion against a stub which exits successfully and
        writes the given bytes, or nothing at all if output is None
        '''
        from pymol import exporting

        completed = subprocess.CompletedProcess([], 0, stdout='', stderr='')

        def fake_run(args, **kwargs):
            if output is not None:
                with open(args[-1], 'wb') as handle:
                    handle.write(output)
            return completed

        with unittest.mock.patch.object(
                exporting, '_find_usdcat', lambda: 'usdcat'), \
                unittest.mock.patch.object(
                    subprocess, 'run', side_effect=fake_run):
            return exporting._convert_usda_to_usdc('#usda 1.0\n')

    def testUSDZConversionOutput(self):
        # a zero exit status alone must not be trusted
        for output, expected in [
                (None, 'no output file'),
                (b'', 'empty'),
                (b'#usda 1.0\n', 'binary USDC')]:
            with self.assertRaises(pymol.CmdException) as caught:
                self.convert_with_fake_usdcat(output)
            self.assertIn(expected, str(caught.exception))

        usdc = b'PXR-USDC' + b'\0' * 8
        self.assertEqual(self.convert_with_fake_usdcat(usdc), usdc)

    def testUSDZPackageAlignment(self):
        from pymol import exporting

        contents = b'PXR-USDC' + b'\0' * 120

        with testing.mktemp('.usdz') as filename:
            exporting._write_usdz_package(filename, 'scene.usdc', contents)

            with zipfile.ZipFile(filename) as archive:
                info = archive.getinfo('scene.usdc')
                self.assertEqual(archive.read('scene.usdc'), contents)

            with open(filename, 'rb') as handle:
                handle.seek(info.header_offset)
                header = struct.unpack('<IHHHHHIIIHH', handle.read(30))

            data_offset = info.header_offset + 30 + header[-2] + header[-1]
            self.assertEqual(data_offset % 64, 0)

    def testUSDZPackageTooLarge(self):
        from pymol import exporting

        # a ZIP64 local header would carry a second extra field and shift the
        # payload off the 64-byte boundary
        with testing.mktemp('.usdz') as filename:
            with unittest.mock.patch.object(zipfile, 'ZIP64_LIMIT', 16):
                with self.assertRaises(pymol.CmdException) as caught:
                    exporting._write_usdz_package(
                        filename, 'scene.usdc', b'x' * 64)

        self.assertIn('too large', str(caught.exception))

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
