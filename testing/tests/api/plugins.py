'''
Plugin contract without a GUI: discovery in a startup directory, import as
pymol.plugins.startup.<name>, __init_plugin__ with the headless app object
(menu calls are no-ops), and cmd.extend registration tracking.
'''

import os
import shutil
import tempfile

from pymol import cmd, plugins, testing

PLUGIN = '''
from pymol import cmd, plugins

def __init_plugin__(app=None):
    plugins.addmenuitem('Test Plugin|Hello', lambda: None)
    cmd.extend('_pymolish_testplugin_hello', lambda: 'hello')
'''


class TestPlugins(testing.PyMOLTestCase):

    def test_headless_plugin_load(self):
        tmp = tempfile.mkdtemp()
        orig = plugins.get_startup_path(True)
        try:
            with open(os.path.join(tmp, 'pymolish_testplugin.py'), 'w') as f:
                f.write(PLUGIN)
            plugins.set_startup_path([tmp], autosave=False)

            found = plugins.findPlugins([tmp])
            self.assertIn('pymolish_testplugin', found)

            info = plugins.PluginInfo('pymolish_testplugin',
                                      found['pymolish_testplugin'],
                                      'pymol.plugins.startup.pymolish_testplugin')
            self.assertTrue(info.load(plugins.get_pmgapp()))
            self.assertTrue(info.loaded)
            self.assertIn('_pymolish_testplugin_hello', info.commands)
            self.assertEqual(cmd.keyword['_pymolish_testplugin_hello'][0](), 'hello')
        finally:
            plugins.set_startup_path(orig, autosave=False)
            shutil.rmtree(tmp)
