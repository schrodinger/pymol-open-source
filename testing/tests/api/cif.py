from pymol import cmd, testing, stored

@testing.requires_version('1.7.1')
class TestCIF(testing.PyMOLTestCase):

    def testPYMOL1737(self):
        cmd.load(self.datafile('ice_IV.cif'))
        self.assertTrue('ice_IV' in cmd.get_object_list())

    def testPYMOL1533(self):
        cmd.load(self.datafile('1v5a-3models.cif'))
        self.assertEqual(cmd.count_states(), 3)

    @testing.requires_version('1.7.7')
    def test_models_diff_atom_count(self):
        cmd.load(self.datafile('1v5a-2models-1st-trunc.cif'), 'm1')
        self.assertEqual(cmd.count_states(), 2)
        self.assertEqual(cmd.count_atoms(), 387)
        self.assertEqual(cmd.count_atoms('state 1'), 139)

    @testing.requires_version('1.7.7')
    def test_chemical_conn_bond(self):
        cmd.load(self.datafile('1519159.cif'), 'm1')
        model = cmd.get_model()
        self.assertEqual(len(model.atom), 26)
        self.assertEqual(len(model.bond), 35)

    @testing.requires_version('1.7.7')
    def test_CCDC(self):
        cmd.load(self.datafile('CAFINE.cif'), 'm1')
        self.assertEqual(cmd.count_atoms(), 25)
        sym = cmd.get_symmetry('m1')
        self.assertEqual(sym[6], 'P 21/a')
        self.assertArrayEqual(sym[:6], [14.8, 16.7, 3.97, 90., 97., 90.], delta=1e-4)

    @testing.requires_version('1.7.7')
    def test_components_multiplex(self):
        cmd.load(self.datafile('components-000-001-002.cif'))
        for name, natoms, nbonds in [
                ('002', 67, 67),
                ('001', 87, 90),
                ('000', 9, 8),
                ]:
            model = cmd.get_model(name)
            self.assertEqual(len(model.atom), natoms)
            self.assertEqual(len(model.bond), nbonds)
        # for 000, test bond orders
        self.assertEqual([b.order for b in model.bond], [1, 2, 1, 1, 1, 1, 1, 1])
