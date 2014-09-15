from pymol import cmd, testing, stored, callback

@testing.requires_version('1.7.3.0')
class TestCallback(testing.PyMOLTestCase):

    def testPseSupport(self):
        cmd.load_callback(callback.Callback(), 'c1')
        with testing.mktemp('tmp.pse') as session_filename:
            cmd.save(session_filename)
            cmd.delete('*')
            cmd.load(session_filename)
            self.assertTrue('c1' in cmd.get_names())
