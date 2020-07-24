import sys
import unittest
from pymol import cmd, testing, stored

try:
    from io import StringIO
    from unittest.mock import patch
    mock_not_available = False
except ImportError:
    mock_not_available = True


def func_with_indented_help():
    '''
    USAGE

        foo

    SEE ALSO

        https://github.com/schrodinger/pymol-open-source/issues/116
    '''


cmd.extend('func_with_indented_help', func_with_indented_help)


@unittest.skipIf(mock_not_available, "unittest.mock not available")
class TestHelping(testing.PyMOLTestCase):
    def testApi(self):
        with patch('sys.stdout', new=StringIO()) as out:
            cmd.api("color")
            self.assertTrue('API: pymol.viewing.color' in out.getvalue())

    def testHelp(self):
        with patch('sys.stdout', new=StringIO()) as out:
            cmd.help('color')
            self.assertTrue('USAGE\n\n    color color' in out.getvalue())

    @testing.requires_version('2.5')
    def testHelp_dedent(self):
        with patch('sys.stdout', new=StringIO()) as out:
            cmd.help('func_with_indented_help')
            self.assertTrue('USAGE\n\n    foo\n\nSEE' in out.getvalue())

    @testing.requires_version('2.4')
    @testing.requires('incentive')
    def testHelpSetting(self):
        out = cmd.help_setting('transparency')
        self.assertTrue('controls surface transparency' in out)
