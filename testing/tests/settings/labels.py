'''
Label settings
'''

from pymol import cmd, testing, stored, invocation

def getQuadrantSlice(img, x, y):
    '''
    Get indices to slice out a quadrant from an array.

    @param x: -1 or 1
    @param y: -1 or 1

     +-----+-----+
     |-1, 1| 1, 1|
     +-----+-----+
     |-1,-1| 1,-1|
     +-----+-----+

    '''
    h, w = img.shape[:2]
    w2, h2 = w // 2, h // 2
    x0 = 0 if x < 0 else w2
    y0 = 0 if y > 0 else h2
    return slice(y0, y0 + h2), slice(x0, x0 + w2)

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

    def assertColorInQuadrant(self, x, y, w, h, color='white'):
        '''
        Assert that the given color is found in the given (x,y) quadrant
        but not in the other three quadrants.
        '''
        img = self.get_imagearray(width=w, height=h)
        idx = getQuadrantSlice(img, x, y)

        self.assertImageHasColor(color, img[idx],
                delta=1,  # observed 254 instead of 255 for white
                msg='no %s in Q(%d,%d)' % (color, x, y))

        img[idx] = 0

        self.assertImageHasNotColor(color, img,
                msg='%s outside of Q(%d,%d)' % (color, x, y))

    @testing.requires('incentive')
    @testing.requires_version('1.8.7')
    def testLabelRelativeMode(self):
        # +-----+-----+
        # |     |     |
        # +-----+-----+
        # | XYZ |     |
        # +-----+-----|

        self.ambientOnly()
        cmd.set('opaque_background')
        cmd.set('label_font_id', 7)  # bold
        cmd.set('label_color', 'white')
        cmd.set('label_size', -1)  # Angstrom sized
        cmd.pseudoatom('m1', label="XYZ")
        cmd.zoom()

        width = 400
        height = 200
        w4 = width // 4
        h4 = height // 4

        cmd.viewport(width, height)

        for x in (-1, 1):
            for y in (-1, 1):
                # relative position
                cmd.set('label_relative_mode', 1)
                cmd.set('label_screen_point', (0.5 * x, 0.5 * y, 0))
                self.assertColorInQuadrant(x, y, width, height)

                # absolute position in pixels
                cmd.set('label_relative_mode', 2)
                cmd.set('label_screen_point', ((2 + x) * w4, (2 + y) * h4, 0))
                self.assertColorInQuadrant(x, y, width, height)
