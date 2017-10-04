'''
PyMOLWin.exe launch tests
'''

import os
import time
import subprocess
import sys
import unittest
from distutils.spawn import find_executable
from pymol import cmd, CmdException, testing, stored

@unittest.skipIf(not sys.platform.startswith('win'), 'windows only')
class TestPyMOLWin(testing.PyMOLTestCase):

    def _get_pymolwin(self):
        pymolwin = os.path.join(os.path.dirname(sys.executable), 'PyMOLWin.exe')

        if os.path.exists(pymolwin):
            return pymolwin

        self.skipTest('missing ' + pymolwin)

    @testing.requires_version('2.0')
    def test_wchar(self):
        '''
        Wide-char file names
        '''
        pymolwin = self._get_pymolwin()

        with testing.mkdtemp() as base:
            dirname = u'Japanese \u65e5\u672c\u8a9e'
            dirname = os.path.join(base, dirname)
            os.mkdir(dirname)

            out_filename = base + '\\out.pdb'
            bat_filename = base + '\\in.bat'
            pml_filename = dirname + u'\\in.pml'
            in_filename = dirname + u'\\in.pdb'
            command = u'"%s" -kcq "%s" "%s"' % (pymolwin, in_filename, pml_filename)

            cmd.fragment('gly')
            cmd.save(in_filename)
            cmd.delete('*')

            with open(pml_filename, 'w') as handle:
                handle.write('save %s\n' % out_filename)

            with open(bat_filename, 'wb') as handle:
                handle.write(b'chcp 65001\r\n')  # UTF-8 code page
                handle.write(command.encode('utf-8') + b'\r\n')
                handle.write(b'chcp 437\r\n')

            subprocess.call([bat_filename])

            # wait for spawned process
            time.sleep(0.5)

            cmd.load(out_filename)
            self.assertTrue(7, cmd.count_atoms())

    def test_unc(self):
        '''
        UNC-Paths (PYMOL-2954)
        '''
        pymolwin = self._get_pymolwin()

        with testing.mkdtemp() as base:
            if base[1] != ':':
                self.skipTest('no drive letter: ' + base)
                return

            if base.startswith('\\\\'):
                uncpath = base
            else:
                uncpath = '\\\\localhost\\%s$%s' % (base[0], base[2:])

            pml_filename = uncpath + '\\in.pml'
            out_filename = uncpath + '\\out.pdb'

            with open(pml_filename, 'w') as handle:
                handle.write('fragment gly\n')
                handle.write('save %s\n' % out_filename)

            subprocess.call([pymolwin, '+2', '-kcq', pml_filename])

            # wait for spawned process
            time.sleep(0.5)

            cmd.load(out_filename)
            self.assertTrue(7, cmd.count_atoms())
