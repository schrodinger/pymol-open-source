from pymol import cmd, testing, stored, CmdException

class TestQuerying(testing.PyMOLTestCase):

    def _testMeasure(self, func, getfunc, natom, value):
        name = "m_obj"
        cmd.fragment("gly")
        
        r = func(name, *["index %d" % (i+1) for i in range(natom)])
        self.assertAlmostEqual(r, value, delta=1e-2)
        self.assertTrue(name in cmd.get_names())
        
        r = getfunc(*["index %d" % (i+1) for i in range(natom)])
        self.assertAlmostEqual(r, value, delta=1e-2)

    def testAngle(self):
        self._testMeasure(cmd.angle, cmd.get_angle, 3, 110.99)

    def testAutoMeasure(self):
        cmd.fragment("gly")

        for (N, prefix) in enumerate(["dist", "angle", "dihedral"], 2):
            cmd.set("dist_counter", 0)
            cmd.edit(*["index %d" % (i + 1) for i in range(N)])
            cmd.auto_measure()
            self.assertTrue(prefix + "01" in cmd.get_names(), prefix)

    def testCountAtoms(self):
        # used in many tests
        pass

    def testCountFrames(self):
        self.assertEqual(cmd.count_frames(), 0)
        cmd.pseudoatom("m1")
        self.assertEqual(cmd.count_frames(), 1)
        cmd.mset("1x3")
        self.assertEqual(cmd.count_frames(), 3)

    def testCountStates(self):
        self.assertEqual(cmd.count_states(), 0)
        cmd.pseudoatom("m1")
        self.assertEqual(cmd.count_states(), 1)
        cmd.create("m1", "m1", 1, 2)
        self.assertEqual(cmd.count_states(), 2)

    def testDihedral(self):
        self._testMeasure(cmd.dihedral, cmd.get_dihedral, 4, 0.01318)

    def testDistance(self):
        self._testMeasure(cmd.distance, cmd.get_distance, 2, 1.46011)

    def testFindPairs(self):
        # mode 0
        cmd.fragment("gly")
        pairs = cmd.find_pairs("hydro", "elem O")
        pairs_ref = [
            (('gly', 5), ('gly', 4)),
            (('gly', 6), ('gly', 4)),
            (('gly', 7), ('gly', 4)),
        ]
        self.assertEqual(sorted(pairs), pairs_ref)

        # mode 1
        cmd.fab("GGGGG", "m2", ss=1)
        pairs = cmd.find_pairs("m2 & donor", "m2 & acceptor", mode=1)
        self.assertEqual(pairs, [(('m2', 29), ('m2', 4))])

    def testGetAngle(self):
        # see testAngle
        pass

    def testGetArea(self):
        cmd.fragment("gly")
        r = cmd.get_area(load_b=1)
        b_list = []
        cmd.iterate("elem O", "b_list.append(b)", space=locals())
        self.assertAlmostEqual(r, 82.505165, delta=1e-2)
        self.assertAlmostEqual(b_list[0], 15.47754, delta=1e-2)
        cmd.set("dot_solvent")
        self.assertAlmostEqual(cmd.get_area(), 200.145, delta=1e-2)
        cmd.set("dot_density", 4)
        self.assertAlmostEqual(cmd.get_area(), 200.888, delta=1e-2)

    def testGetAtomCoords(self):
        cmd.fragment("gly")
        coords = cmd.get_atom_coords("elem O")
        coords_ref = [0.545436, -0.974895, 1.498943]
        self.assertArrayEqual(coords, coords_ref, delta=1e-4)
        self.assertRaises(CmdException, cmd.get_atom_coords, "*")

    def testGetChains(self):
        cmd.pseudoatom("m1", chain="C")
        self.assertEqual(cmd.get_chains(), ["C"])

    def testGetColorIndex(self):
        self.assertEqual(cmd.get_color_index("blue"), 2)

    def testGetColorIndices(self):
        r = cmd.get_color_indices()
        self.assertTrue(set(r).issuperset([('white', 0), ('black', 1), ('blue', 2)]))

    def testGetColorTuple(self):
        self.assertEqual((0.0, 0.0, 1.0), cmd.get_color_tuple("blue"))

    def testGetDihedral(self):
        # see testDihedral
        pass

    def testGetDistance(self):
        # see testDistance
        pass

    def testGetDragObjectName(self):
        cmd.get_drag_object_name
        self.skipTest("TODO")

    def testGetExtent(self):
        cmd.get_extent
        self.skipTest("TODO")

    def testGetIdtf(self):
        cmd.get_idtf
        self.skipTest("TODO")

    def testGetLegalName(self):
        cmd.get_legal_name
        self.skipTest("TODO")

    def testGetModalDraw(self):
        cmd.get_modal_draw
        self.skipTest("TODO")

    def testGetModel(self):
        cmd.get_model
        self.skipTest("TODO")

    def testGetMovieLength(self):
        cmd.get_movie_length
        self.skipTest("TODO")

    def testGetMovieLocked(self):
        cmd.get_movie_locked
        self.skipTest("TODO")

    def testGetMtlObj(self):
        cmd.get_mtl_obj
        self.skipTest("TODO")

    def testGetNames(self):
        cmd.get_names
        self.skipTest("TODO")

    def testGetNamesOfType(self):
        cmd.get_names_of_type
        self.skipTest("TODO")

    def testGetObjectColorIndex(self):
        cmd.get_object_color_index
        self.skipTest("TODO")

    def testGetObjectList(self):
        cmd.get_object_list
        self.skipTest("TODO")

    def testGetObjectMatrix(self):
        cmd.get_object_matrix
        self.skipTest("TODO")

    def testGetPhipsi(self):
        cmd.get_phipsi
        self.skipTest("TODO")

    def testGetPosition(self):
        cmd.get_position
        self.skipTest("TODO")

    def testGetPovray(self):
        cmd.get_povray
        self.skipTest("TODO")

    def testGetRawAlignment(self):
        cmd.get_raw_alignment
        self.skipTest("TODO")

    def testGetRenderer(self):
        cmd.get_renderer
        self.skipTest("TODO")

    def testGetSymmetry(self):
        cmd.get_symmetry
        self.skipTest("TODO")

    def testGetTitle(self):
        cmd.get_title
        self.skipTest("TODO")

    def testGetType(self):
        cmd.get_type
        self.skipTest("TODO")

    def testGetUnusedName(self):
        cmd.get_unused_name
        self.skipTest("TODO")

    def testGetVersion(self):
        cmd.get_version
        self.skipTest("TODO")

    def testGetVolumeField(self):
        cmd.get_volume_field
        self.skipTest("TODO")

    def testGetVolumeHistogram(self):
        cmd.get_volume_histogram
        self.skipTest("TODO")

    def testGetVolumeIsUpdated(self):
        cmd.get_volume_is_updated
        self.skipTest("TODO")

    def testGetVolumeRamp(self):
        cmd.get_volume_ramp
        self.skipTest("TODO")

    def testGetVrml(self):
        cmd.get_vrml
        self.skipTest("TODO")

    def testIdAtom(self):
        cmd.fragment('gly', 'm1')
        self.assertEqual(cmd.id_atom('ID 3'), 3)
        self.assertEqual(cmd.id_atom('ID 3', 1), ('m1', 3))
        cmd.feedback("disable", "cmd", "errors")
        self.assertRaises(CmdException, cmd.id_atom, 'ID 3+4')
        self.assertRaises(CmdException, cmd.id_atom, 'ID 100')

    def testIdentify(self):
        cmd.fragment('gly', 'm1')
        cmd.select('s1', 'ID 3+4')
        r = cmd.identify('s1')
        self.assertItemsEqual(r, [3, 4])
        r = cmd.identify('s1', 1)
        self.assertItemsEqual(r, [('m1', 3), ('m1', 4)])

    def testIndex(self):
        cmd.fragment('gly', 'm1')
        r = cmd.index('index 3+4')
        self.assertItemsEqual(r, [('m1', 3), ('m1', 4)])

    def testOverlap(self):
        # unsupported command
        cmd.pseudoatom('m1', vdw=2.0)
        cmd.pseudoatom('m2', vdw=2.0, pos=(2., 0., 0.))
        r = cmd.overlap('m1', 'm2')
        # this is actually half the vdw overlap
        self.assertEqual(r, 1.0)

    def testPhiPsi(self):
        cmd.phi_psi
        self.skipTest("TODO")

    def testSetVolumeRamp(self):
        cmd.set_volume_ramp
        self.skipTest("TODO")

