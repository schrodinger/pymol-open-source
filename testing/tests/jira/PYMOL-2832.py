from pymol import cmd, testing

@testing.requires_version('1.8.4.1')
@testing.requires('incentive')
class TestDiscreteRibbonColors(testing.PyMOLTestCase):

    @testing.foreach(1, 0)
    def test(self, ribbon_as_cylinders):
        cmd.set('ribbon_as_cylinders', ribbon_as_cylinders)
        cmd.set('ribbon_width', 8)
        cmd.fab('AG')
        cmd.color('red', 'resn ALA')
        cmd.color('blue', 'resn GLY')
        cmd.show_as('ribbon')
        cmd.orient()

        self.ambientOnly()
        cmd.viewport(100, 100)
        cmd.draw(100, 100)
        img = self.get_imagearray()
        self.assertImageHasColor('red', img)
        self.assertImageHasColor('blue', img)
