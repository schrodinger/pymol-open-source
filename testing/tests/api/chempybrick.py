from chempy import cpv
from pymol import cmd, testing, stored

class TestChempyBrick(testing.PyMOLTestCase):

    def test(self):
        try:
            import numpy
        except ImportError:
            self.skip('requires numpy')

        from chempy.brick import Brick

        spacing = (.5, .5, .5)
        origin = (3., 4., 5.)
        shape = (10, 10, 10)
        range_ =  [a * (b - 1) for (a, b) in zip (spacing, shape)]

        data = numpy.zeros(shape, float)

        brick = Brick.from_numpy(data, spacing, origin)
        self.assertArrayEqual(brick.range, range_)

        cmd.load_brick(brick, 'map')

        extent = cmd.get_extent('map')
        self.assertArrayEqual(extent, [origin, cpv.add(origin, range_)])
