python
import sys, os
sys.argv = ['pymol', '--run', 'all', '--offline']
cmd.do('run ' + os.path.join(os.path.dirname(__script__), 'testing.py'))
python end
