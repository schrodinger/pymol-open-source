from pymol import testing
from pymol import constants

class TestConstants(testing.PyMOLTestCase):
    def test_safe_alpha_list_eval(self):
        f = constants.safe_alpha_list_eval
        # expected use case
        self.assertEqual(f('[blue, red, green]'), ["blue", "red", "green"])
        # single item
        self.assertEqual(f('[AB]'), ["AB"])
        # trailing comma
        self.assertEqual(f('[A,]'), ["A"])
        # tuple
        self.assertEqual(f("('A',)"), ('A',))
        # single value
        self.assertEqual(f("('A')"), 'A')
        self.assertEqual(f("'A'"), 'A')
        # with quotes
        self.assertEqual(f('["AB", \'cd\']'), ["AB", "cd"])
        # str, int, non-alnum characters
        self.assertEqual(f('[A, 1, 2+3, B:C, D E, "F G"]'),
                         ["A", 1, 23, "BC", "DE", "FG"])

    @testing.requires_version('1.6')
    def test_safe_eval(self):
        f = constants.safe_eval
        self.assertEqual(f('[AB]'), ["AB"])
        self.assertEqual(f('[A,]'), ["A"])
        self.assertEqual(f("('A',)"), ('A',))
        self.assertEqual(f("('A')"), 'A')
        self.assertEqual(f("'A'"), 'A')
        self.assertEqual(f('["AB", \'cd\']'), ["AB", "cd"])
        self.assertEqual(f('[A, 1, 2+3, "F G"]'), ["A", 1, 5, "F G"])
        self.assertEqual(f('()'), ())
        self.assertEqual(f('[]'), [])
        self.assertEqual(f('{}'), {})
        self.assertEqual(f('{foo: bar}'), {"foo": "bar"})
        self.assertEqual(f('1 + 2'), 3)
        with self.assertRaises(TypeError):
            f('__import__("os").times()')
