import os, sys
from pymol import cmd, testing, stored, invocation

class TestInvocation(testing.PyMOLTestCase):

    @testing.foreach.zip(('.pymolrc.py', '.pymolrc', 'pymolrc.pml', 'pymolrcextra.py'))
    def test_get_user_config(self, basename):
        '''
        pymolrc finding test
        '''

        orig = dict((key, os.getenv(key)) for key in ['HOME', 'HOMEDRIVE', 'HOMEPATH', 'PYMOL_PATH'])

        with testing.mkdtemp() as tmpdir:
            os.chdir(tmpdir)
            os.mkdir('abc')

            pwd = os.getcwd()

            if sys.platform.startswith('win'):
                if pwd[1] != ':':
                    raise ValueError('no drive letter')
                os.environ.pop('HOME', None)
                os.environ['HOMEDRIVE'] = pwd[:2]
                os.environ['HOMEPATH'] = pwd[2:] + r'\abc'
            else:
                os.environ['HOME'] = pwd + '/abc'

            pymolrc = os.path.join(pwd, 'abc', basename) 
            with open(pymolrc, 'w') as handle:
                print >> handle, '# hello'

            a = invocation.get_user_config()
            self.assertEqual(pymolrc, os.path.normpath(a[0]))

            pymolrc = os.path.join(pwd, basename)
            with open(pymolrc, 'w') as handle:
                print >> handle, '# hello'

            a = invocation.get_user_config()
            self.assertEqual(pymolrc, os.path.normpath(a[0]))

        for key, value in orig.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

