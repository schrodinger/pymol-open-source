'''
CCP4/MRC ORIGIN
'''

from chempy import cpv
from pymol import cmd, testing

EDGELENGTHS = (5.85, 5.85, 5.85)
ORIGINS = {
    # only SKWTRN
    "h2o-elf-skwtrn.ccp4": (1.2, 3.4, 5.6),
    # NSTART and SKWTRN, additive
    "h2o-elf-nstart-skwtrn.ccp4": (1.8, 3.85, 6.35),
    # NSTART and ORIGIN, not additive (NSTART is ignored)
    "h2o-elf-nstart-origin.mrc": (1.2, 3.4, 5.6),
    # only NSTART (not ignored)
    "h2o-elf-nstart.ccp4": (0.6, 0.45, 0.75),
}

@testing.requires_version('2.4')
class Test3260(testing.PyMOLTestCase):

    @testing.foreach(*ORIGINS)
    def test(self, filename):
        cmd.load(self.datafile(filename), 'map')

        ori = ORIGINS[filename]
        ext = [ori, cpv.add(ori, EDGELENGTHS)]

        self.assertArrayEqual(cmd.get_extent('map'), ext, delta=1e-3)

        with testing.mktemp('.map') as filename:
            cmd.save(filename, 'map')
            cmd.delete('*')
            cmd.load(filename, 'map')

        self.assertArrayEqual(cmd.get_extent('map'), ext, delta=1e-3)
