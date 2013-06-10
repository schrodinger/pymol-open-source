from pymol import cmd, testing, stored

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
        cmd.load_model
        self.skipTest("TODO")

    def testLoadObject(self):
        cmd.load_object
        self.skipTest("TODO")

    def testLoadRaw(self):
        cmd.load_raw
        self.skipTest("TODO")

    def testLoadTraj(self):
        cmd.load_traj
        self.skipTest("TODO")

    def testReadMmodstr(self):
        cmd.read_mmodstr
        self.skipTest("TODO")

    def testReadMolstr(self):
        cmd.read_molstr
        self.skipTest("TODO")

    def testReadPdbstr(self):
        cmd.read_pdbstr
        self.skipTest("TODO")

    def testReadSdfstr(self):
        cmd.read_sdfstr
        self.skipTest("TODO")

    def testReadXplorstr(self):
        cmd.read_xplorstr
        self.skipTest("TODO")

    def testSetSession(self):
        cmd.set_session
        self.skipTest("TODO")

    def testSpace(self):
        cmd.space
        self.skipTest("TODO")

