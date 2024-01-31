'''
Non-ASCII paths
'''

import os
from pymol import cmd, CmdException, testing, stored

@testing.requires_version('1.8.7')
class TestUnicodePaths(testing.PyMOLTestCase):
    def test(self):
        cmd.fragment('gly', 'm1')

        with testing.mkdtemp() as base:
            for dirname in [
                    u'Japanese \u65e5\u672c\u8a9e',
                    u'Umlauts \xc4\xd6\xdc',
                ]:
                dirname = os.path.join(base, dirname)
                os.mkdir(dirname)

                # PNG image
                filename = os.path.join(dirname, u'image.png')
                cmd.png(filename, 10, 10, ray=1)
                self.assertTrue(os.path.isfile(filename))

                # Data files
                for filename in [
                        u'm1.pdb',
                        u'session.pse',
                    ]:
                    filename = os.path.join(dirname, filename)

                    cmd.save(filename)
                    self.assertTrue(os.path.isfile(filename))

                    cmd.delete('*')
                    cmd.load(filename)
                    self.assertEqual(cmd.get_object_list(), ['m1'])
