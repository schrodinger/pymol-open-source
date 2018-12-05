from pymol import cmd, stored, testing

@testing.requires_version('2.3')
class TestParserNesting(testing.PyMOLTestCase):
    def test(self):
        # basic
        stored.x = []
        cmd.do(
            'cmd.do("'
            '  stored.x.append(1);'
            '  stored.x.append(2)'
            '");'
            'stored.x.append(3)',
            echo=0)
        self.assertEqual(stored.x, [1, 2, 3])

        # extra level of nesting
        stored.x = []
        cmd.do(
            'cmd.do("'
            '  stored.x.append(1);'
            "  cmd.do('"
            '    stored.x.append(2);'
            '    stored.x.append(3)'
            "  ');"
            '  stored.x.append(4)'
            '");'
            'stored.x.append(5)',
            echo=0)
        self.assertEqual(stored.x, [1, 2, 3, 4, 5])
