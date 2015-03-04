'''
test stick rep settings
'''

from pymol import cmd, testing

class TestStick(testing.PyMOLTestCase):

    def testStickBall(self):
        cmd.viewport(100, 100)

        self.ambientOnly()

        cmd.set('stick_ball')
        cmd.set('stick_ball_ratio', 2.0)
        cmd.set('stick_ball_color', 'blue')
        cmd.set('stick_color', 'red')

        cmd.fragment('gly')
        cmd.orient()

        cmd.color('green')
        cmd.show_as('sticks')

        img = self.get_imagearray()
        self.assertImageHasColor('blue', img)
        self.assertImageHasColor('red', img)
        self.assertImageHasNotColor('green', img)
