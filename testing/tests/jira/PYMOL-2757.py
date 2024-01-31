'''
crash with mesh w/o map
'''

from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.8.0.6') # 1.8.0.7
class TestPYMOL2757(testing.PyMOLTestCase):

    def test(self):
        cmd.viewport(100, 100)

        # make map
        cmd.fragment('gly', 'm1')
        cmd.set('gaussian_b_floor', 30)
        cmd.set('mesh_width', 5)
        cmd.map_new('map')
        cmd.delete('m1')

        # make mesh 
        cmd.isomesh('mesh', 'map')

        # check mesh presence by color
        meshcolor = 'red'
        cmd.color(meshcolor, 'mesh')
        self.ambientOnly()
        self.assertImageHasColor(meshcolor)

        # continue without map
        cmd.delete('map')

        with testing.mktemp('.pse') as filename:
            cmd.save(filename)

            cmd.delete('*')
            self.assertImageHasNotColor(meshcolor)

            cmd.load(filename)
            self.assertImageHasColor(meshcolor)

