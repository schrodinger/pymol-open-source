from pymol import cmd, testing, stored

pdbstr = '''ATOM      1  N   GLY    22      -1.195   0.201  -0.206  1.00  0.00           N  
ATOM      2  CA  GLY    22       0.230   0.318  -0.502  1.00  0.00           C  
ATOM      3  C   GLY    22       1.059  -0.390   0.542  1.00  0.00           C  
ATOM      4  O   GLY    22       0.545  -0.975   1.499  1.00  0.00           O  
ATOM      5  H   GLY    22      -1.558  -0.333   0.660  1.00  0.00           H  
ATOM      6 3HA  GLY    22       0.482   1.337  -0.514  0.00  0.00           H  
ATOM      7  HA  GLY    22       0.434  -0.159  -1.479  1.00  0.00           H  
END
'''

molstr = '''glycine
  ChemPy            3D                             0

  7  6  0  0  1  0  0  0  0  0999 V2000
   -1.1946    0.2011   -0.2061 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2304    0.3181   -0.5021 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0594   -0.3899    0.5419 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5454   -0.9749    1.4989 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5576   -0.3329    0.6599 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4824    1.3374   -0.5137 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4344   -0.1589   -1.4791 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  3  4  2  0  0  0  0
M  END
'''

mmodstr = '''     7  gly
  34     2 1     5 1     0 0     0 0     0 0     0 0   -1.195000    0.201000   -0.206000    22G   38  0.00000  0.00000 GLY   N
   3     1 1     3 1     6 1     7 1     0 0     0 0    0.230000    0.318000   -0.502000    22G    2  0.00000  0.00000 GLY   CA
   7     2 1     4 2     0 0     0 0     0 0     0 0    1.059000   -0.390000    0.542000    22G    2  0.00000  0.00000 GLY   C
  15     3 2     0 0     0 0     0 0     0 0     0 0    0.545000   -0.975000    1.499000    22G   75  0.00000  0.00000 GLY   O
  44     1 1     0 0     0 0     0 0     0 0     0 0   -1.558000   -0.333000    0.660000    22G   21  0.00000  0.00000 GLY   H
  41     2 1     0 0     0 0     0 0     0 0     0 0    0.482000    1.337000   -0.514000    22G   21  0.00000  0.00000 GLY  3HA
  41     2 1     0 0     0 0     0 0     0 0     0 0    0.434000   -0.159000   -1.479000    22G   21  0.00000  0.00000 GLY   HA
'''

class TestImporting(testing.PyMOLTestCase):

    @testing.requires('network', 'no_run_all')
    def testFetch(self):
        with testing.mkdtemp() as fetch_path:
            names = []
            cmd.set('fetch_path', fetch_path)
            cmd.set('fetch_host', 'pdb pdbe')

            cmd.fetch('1avy', '1avy1', type='pdb1')
            names += ['1avy1']
            self.assertItemsEqual(cmd.get_names(), names)
            self.assertEqual(cmd.count_states('1avy1'), 3)

            cmd.fetch('1avy', type='2fofc')
            names += ['1avy_2fofc']
            self.assertItemsEqual(cmd.get_names(), names)
            self.assertEqual(cmd.get_type('1avy_2fofc'), 'object:map')

    def testFetchLocal(self):
        import urlparse
        with testing.mkdtemp() as fetch_path:
            names = []
            cmd.set('fetch_path', fetch_path)
            cmd.set('fetch_host', urlparse.urlunsplit(['file', '',
                self.datafile('pdb.mirror'), '', '']))

            cmd.fetch('1avy')
            names += ['1avy']
            self.assertItemsEqual(cmd.get_names(), names)

            cmd.fetch('1avyB')
            names += ['1avyB']
            self.assertItemsEqual(cmd.get_names(), names)
            self.assertEqual(cmd.get_chains('1avyB'), ['B'])

            cmd.fetch('1aq5', multiplex=1)
            names += ['1aq5_%04d' % (i+1) for i in range(20)]
            self.assertItemsEqual(cmd.get_names(), names)

    def testFinishObject(self):
        cmd.finish_object
        self.skipTest("TODO")

    def testLoad(self):
        cmd.load
        self.skipTest("TODO")

    def testLoadBrick(self):
        cmd.load_brick
        self.skipTest("TODO")

    def testLoadCallback(self):
        cmd.load_callback
        self.skipTest("TODO")

    def testLoadCgo(self):
        cmd.load_cgo
        self.skipTest("TODO")

    def testLoadEmbedded(self):
        cmd.load_embedded
        self.skipTest("TODO")

    def testLoadMap(self):
        cmd.load_map
        self.skipTest("TODO")

    def testLoadModel(self):
        cmd.fragment('gly')
        m = cmd.get_model()
        cmd.delete('*')
        cmd.load_model(m, 'm1')
        self.assertEqual(7, cmd.count_atoms())

    def testLoadObject(self):
        cmd.load_object
        self.skipTest("TODO")

    def testLoadRaw(self):
        cmd.load_raw
        self.skipTest("TODO")

    def testLoadTraj(self):
        pdbfile = self.datafile("sampletrajectory.pdb")
        dcdfile = self.datafile("sampletrajectory.dcd")

        cmd.load(pdbfile)
        cmd.load_traj(dcdfile)
        self.assertEqual(501, cmd.count_states())
        cmd.delete('*')

        cmd.load(pdbfile)
        cmd.load_traj(dcdfile, state=1, interval=5, max=20)
        self.assertEqual(20, cmd.count_states())
        cmd.delete('*')

        cmd.load(pdbfile)
        cmd.load_traj(dcdfile, state=1, start=31, stop=40)
        self.assertEqual(10, cmd.count_states())
        cmd.delete('*')

        cmd.load(pdbfile)
        cmd.load_traj(dcdfile, state=1, stop=30, average=3)
        self.assertEqual(10, cmd.count_states())
        cmd.delete('*')

    def testReadMmodstr(self):
        cmd.read_mmodstr(mmodstr, 'm1')
        self.assertEqual(7, cmd.count_atoms())

    def testReadMolstr(self):
        cmd.read_molstr(molstr, 'm1')
        self.assertEqual(7, cmd.count_atoms())

    def testReadPdbstr(self):
        cmd.read_pdbstr(pdbstr, 'm1')
        self.assertEqual(7, cmd.count_atoms())

        cmd.read_pdbstr(pdbstr, 'm1')
        self.assertEqual(2, cmd.count_states())

    def testReadSdfstr(self):
        sdfstr = molstr + '$$$$\n'

        cmd.read_sdfstr(sdfstr, 'm1')
        self.assertEqual(7, cmd.count_atoms())

    def testReadXplorstr(self):
        cmd.read_xplorstr
        self.skipTest("TODO")

    def testSetSession(self):
        # see TestExporting.testGetSession
        pass

    def testSpace(self):
        cmd.space('pymol', 0.5)
        cmd.space('cmyk', 0.7)
        cmd.space('rgb')

    @testing.requires_version('1.7.3.0')
    def testLoadCoords(self):
        import numpy
        cmd.fragment('gly', 'm1')
        coords = cmd.get_coords('m1')
        coords += 5.0
        cmd.load_coords(coords, 'm1')
        self.assertTrue(numpy.allclose(coords, cmd.get_coords('m1')))

    @testing.requires_version('1.7.3.0')
    def testLoadCoordset(self):
        import numpy
        cmd.fragment('gly', 'm1')
        coords = cmd.get_coordset('m1')
        cmd.load_coordset(coords, 'm1', state=2)
        self.assertEqual(2, cmd.count_states('m1'))

        # data manipulation with copy=0
        cmd.get_coordset('m1', copy=0)[:] += 5.0
        self.assertTrue(numpy.allclose(coords + 5.0, cmd.get_coords('m1')))
