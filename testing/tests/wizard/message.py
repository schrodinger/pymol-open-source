'''
unit tests for pymol.wizard.message
'''

from pymol import cmd, testing

class TestMessage(testing.PyMOLTestCase):

    def test_NonAsciiPrompt(self):
        # github #517: prompt buffer was sized by code points instead of
        # UTF-8 bytes, overflowing the heap for multi-byte characters
        msg = "Unicode test: " + "测试消息" * 2000
        cmd.wizard("message", msg)
        for _ in range(10):
            cmd.refresh_wizard()
        self.assertEqual(cmd.get_wizard().message, [msg])
        cmd.set_wizard()
