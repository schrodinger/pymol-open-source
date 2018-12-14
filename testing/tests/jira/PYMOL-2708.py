'''
crash with bumps and PSE
'''

from pymol import cmd, CmdException, testing, stored

class TestPYMOL2708(testing.PyMOLTestCase):

    # works in (1.7.4 <= version <= 1.7.6) || (version >= 1.8.0.5)
    # broken in (1.8.0.0 <= version <= 1.8.0.4)
    @testing.requires_version('1.8.0.5')
    def testBumpsPSE(self):
        self.ambientOnly()
        cmd.viewport(100, 100)

        # expect some clashes in this helix
        cmd.fab('AAAAA', 'm1', ss=1)
        cmd.orient()

        cmd.set('sculpt_vdw_vis_mode', 1)
        cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
        cmd.sculpt_activate('m1')
        cmd.sculpt_iterate('m1', cycles=0)
        cmd.show_as('cgo')

        img = self.get_imagearray()
        self.assertImageHasColor('0xFF3333', img)

        with testing.mktemp('.pse') as filename:
            cmd.save(filename)
            cmd.load(filename)

        self.assertImageEqual(img)
