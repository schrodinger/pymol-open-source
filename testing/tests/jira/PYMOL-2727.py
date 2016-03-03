'''
missing residues
'''

from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.8.1.0')
class TestPYMOL2727(testing.PyMOLTestCase):

    def _get_seq(self, selection):
        from pymol.exporting import _resn_to_aa as one_letter
        prev_resi = [None]
        seq = []
        def callback(resi, resn):
            if resi == prev_resi[0]:
                return
            prev_resi[0] = resi
            seq.append(one_letter[resn])
        cmd.iterate(selection, 'callback(resi, resn)', space=locals())
        return ''.join(seq)

    @testing.foreach.product(
            (0, 1),
            ('1hbb_pdbx_seq_one_letter_code.cif', '1hbb_entity_poly_seq.cif'),
        )
    def testMissingRes(self, use_auth, filename):
        cmd.set('cif_use_auth', use_auth)
        cmd.load(self.datafile(filename), '1hbb')

        seq = self._get_seq('segi A')
        self.assertEqual(seq, 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR')

        seq = self._get_seq('segi B')
        self.assertEqual(seq, 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH')
        
