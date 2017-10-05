python
import sys, os
_testing_py = os.path.join(os.path.dirname(__script__), 'testing.py')
sys.argv = [_testing_py, '--run', 'all', '--offline', '--no-mmlibs']
cmd.do('run ' + _testing_py)
python end
