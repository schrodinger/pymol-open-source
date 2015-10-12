'''
Crash with names selection of length >= 1024
'''

from pymol import cmd, util, testing, stored

@testing.requires('incentive')
@testing.requires_version('1.7.6.6')
class Test2588(testing.PyMOLTestCase):

    def test(self):
        for i in range(100):
            cmd.pseudoatom('xxxxxxxxxx%04d' % i)

        setting = 'sphere_color'
        color = 'blue'
        name = 'xxxxxxxxxx0005'
        cmd.set(setting, color, name)
        self.assertEqual(cmd.get(setting, name), color)

        # crashes in 1.7.6.0 to 1.7.6.6
        util.color_deep('red', 'all')

        self.assertEqual(cmd.get(setting, name), 'default')
