from pymol.cgo import *
import pymol
from pymol import cmd, testing, stored
import unittest

class TestCGOLighting(testing.PyMOLTestCase):

    @testing.requires('gui')
    @testing.requires_version('1.7.3')
    @testing.foreach((0),(1))
    # if normals are set, then object's cgo_lighting should be set
    def test_cgo_lighting_setting(self, has_normal):
        objname = 'simpletri'
        obj = []
        if has_normal:
            obj.extend([NORMAL, 0., 0., 1.])
        obj.extend([
            BEGIN, TRIANGLES,
            COLOR, 1.0, 0., 0.,
            VERTEX, 0.0, 0.0, 0.0,
            VERTEX, 1.0, 0.0, 0.0,
            VERTEX, 0.0, 1.0, 0.0,
            END
            ])
        cmd.load_cgo(obj,objname,1)
        cgo_lighting = cmd.get('cgo_lighting', objname)
        self.assertEqual(int(cgo_lighting), int(has_normal))

    @testing.requires('gui')
    @testing.requires_version('1.7.3')
    @testing.foreach.product((0,1), (0,1), (0,1), (0,1), (0,1), (0,1))
    def testBackfaceLighting(self, use_shader, has_normal, two_sided_lighting, backface_cull, cgo_lighting, back):
        objname = 'simpletri'
        obj = []
        if has_normal: # specifies a normal
            obj.extend([NORMAL, 0., 0., 1.])
        obj.extend([
            BEGIN, TRIANGLES,
            COLOR, 1.0, 0., 0.,
            VERTEX, 0.0, 0.0, 0.0,
            VERTEX, 1.0, 0.0, 0.0,
            VERTEX, 0.0, 1.0, 0.0,
            END
            ])
        cmd.load_cgo(obj,objname,1)

        cmd.set('use_shader', use_shader)
        cmd.set('two_sided_lighting', two_sided_lighting)
        cmd.set('backface_cull', backface_cull)
        cmd.set('cgo_lighting', cgo_lighting, objname)

        cmd.set('depth_cue', 0)
        cmd.set('antialias', 0)
        cmd.set('light_count', 1)

        if back:
            cmd.turn("y", 180)

        # to see the images
        #        cmd.png("%d__%d_%d_%d_%d_%d_%d_test.png" % (back and backface_cull, use_shader, has_normal, two_sided_lighting, backface_cull, cgo_lighting, back))
        img = self.get_imagearray(width=100, height=100)
        if not two_sided_lighting and not backface_cull and cgo_lighting and back:
            # this is where a shaded back side is rendered
            redmask = (img[:,:,0] > 0)
            self.assertTrue((img[:,:,0][redmask] < 255).any() and
                            (img[:,:,1:3][redmask] == 0).all())
        elif back and backface_cull:
            self.assertImageHasNotColor('red', img=img)
        else:
            self.assertImageHasColor('red', img=img)

