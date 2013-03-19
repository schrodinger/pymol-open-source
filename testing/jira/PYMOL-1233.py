'''
PYMOL-1233 Unable to color isomesh or isodot
'''

import unittest
from pymol import cmd, testing

@testing.requires('gui')
class Test829(testing.PyMOLTestCase):

    def _check_colors(self, *colors):
        for ray, shaders, cylinders in [
                (0,0,0),
                (0,1,0),
                (0,1,1),
                (1,1,1)]:
            cmd.set('use_shaders', shaders)
            cmd.set('mesh_as_cylinders', cylinders)
            img_array = self.get_imagearray(ray=1)
            for color in colors:
                self.assertImageHasColor(color, img_array)

    def test(self):
        cmd.viewport(100,100)

        cmd.set('gaussian_b_floor', 20)
        cmd.set('ambient', 1)
        cmd.set('specular', 0)
        cmd.set('mesh_color', 'blue')
        cmd.set('dot_color', 'green')

        cmd.fragment('gly')
        cmd.map_new('map1')
        cmd.disable('*')

        cmd.isomesh('o1', 'map1')
        cmd.color('red', 'o1')
        cmd.show('cell')
        self._check_colors('red', 'blue')
        cmd.delete('o1')
 
        cmd.isodot('o1', 'map1')
        cmd.color('red', 'o1')
        self._check_colors('green')
        cmd.delete('o1')

        cmd.gradient('o1', 'map1')
        cmd.color('yellow', 'o1')
        self._check_colors('yellow')
        cmd.delete('o1')
