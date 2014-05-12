'''
Label settings
'''

from pymol import cmd, testing, stored, invocation

class TestLabels(testing.PyMOLTestCase):

    @testing.foreach(
        (0, 1),
        (0, 0),
        (1, 0),
    )
    def testLabelPositionZ(self, use_shaders, ray):
        '''
        Test label z position for regular labels
        '''
        if not ray and invocation.options.no_gui:
            self.skipTest('no gui')

        cmd.set('use_shaders', use_shaders)

        cmd.viewport(200, 200)

        cmd.pseudoatom('m1', vdw=10, label='X')
        cmd.zoom(buffer=12)

        cmd.show('spheres')
        cmd.color('blue')
        cmd.set('label_color', 'red')
        cmd.set('label_size', 20)
        cmd.set('label_font_id', 7) # bold
        self.ambientOnly()

        # label outside of sphere
        cmd.set('label_position', [0, 0, 11.1])
        img = self.get_imagearray(ray=ray)
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor('red', img, delta=20)

        # label inside of sphere
        cmd.set('label_position', [0, 0, 10.5])
        img = self.get_imagearray(ray=ray)
        self.assertImageHasColor('blue', img)
        self.assertImageHasNotColor('red', img, delta=20)

