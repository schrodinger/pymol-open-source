'''
unit tests for wizard prompts
'''

from pymol import cmd, testing
from pymol.wizard import Wizard

class TestPrompt(testing.PyMOLTestCase):

    def test_NonAsciiPrompt(self):
        # github #517: prompt buffer was sized by code points instead of
        # UTF-8 bytes, overflowing the heap for multi-byte characters.
        # Plain Wizard instead of the message wizard, which prints the
        # text and fails on non-UTF-8 consoles (e.g. cp1252 on Windows).
        wiz = Wizard()
        wiz.prompt = ["Unicode test: " + "测试消息" * 2000]
        cmd.set_wizard(wiz)
        for _ in range(10):
            cmd.refresh_wizard()
        self.assertIs(cmd.get_wizard(), wiz)
        cmd.set_wizard()
