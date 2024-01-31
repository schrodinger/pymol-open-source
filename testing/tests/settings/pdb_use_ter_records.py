'''
test various settings
'''

import os
from pymol import cmd, testing


@testing.requires_version('2.3')
class TestTER(testing.PyMOLTestCase):

    FILES = [
        'gap.pdb',
        'hetatm-polymer.pdb',
        'multi-model.pdb',
        'no-polymer.pdb',
        'no-solvent.pdb',
        'solvent.pdb',
        'atom-solvent.pdb',
        'two-chains.pdb',
    ]

    @staticmethod
    def ter_indices(lines):
        '''Indices of TER records'''
        return [i for (i, line) in enumerate(lines) if line.startswith('TER ')]

    @testing.foreach.product(FILES, (0, 1))
    def test(self, basename, use_ter):
        cmd.set('pdb_use_ter_records', use_ter)

        filename = self.datafile(os.path.join('ter-records', basename))
        cmd.load(filename, 'm1')

        with open(filename) as handle:
            lines_ref = handle.readlines()

        state = 0 if any(
            line.startswith('MODEL') for line in lines_ref) else -1

        with testing.mktemp('.pdb') as tempfile:
            cmd.save(tempfile, state=state)
            with open(tempfile) as handle:
                lines_out = handle.readlines()

        expected = self.ter_indices(lines_ref) if use_ter else []

        self.assertEqual(expected, self.ter_indices(lines_out))
