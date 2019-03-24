from pymol import cmd, testing

@testing.requires_version('2.2')
@testing.requires('properties')
class TestMaePropertiesIO(testing.PyMOLTestCase):
    def test(self):
        cmd.fragment('gly')

        # atom properties
        cmd.alter('index 2-4', 'p.r_custom_indexhalve = index / 2.0')
        cmd.alter('index 3-5', 'p.neg_index = -index')
        cmd.alter('elem C', 'p["spaced resn name"] = resn + " " + name')

        # object properties
        objprops = [
            ('first', 'Max'),
            ('last', 'Mustermann'),
            ('full name', 'Max Mustermann'),
            ('quoted "name"', 'Max "Mustermann"'),
        ]
        for (key, value) in objprops:
            cmd.set_property(key, value)

        # save/load round-trip
        with testing.mktemp('.mae') as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.set('load_object_props_default', '*')
            cmd.set('load_atom_props_default', '*')
            cmd.load(filename, 'm1')

        # type conversions
        cmd.alter('m1', 'p.neg_index = int(p.s_pymol_neg_index or 0)')

        # atom properties
        self.assertEqual(2, cmd.count_atoms('p.neg_index < -3'))
        self.assertEqual(2, cmd.count_atoms('p.r_custom_indexhalve > 1.0'))

        # changed in 2.3.1: undefined (<>) -> None
        resnname = set()
        cmd.iterate('m1', 'resnname.add(p["s_pymol_spaced resn name"] or "")', space=locals())
        self.assertEqual(resnname, set(['GLY C', 'GLY CA', '']))

        # object properties
        for (key, value) in objprops:
            self.assertEqual(cmd.get_property('s_pymol_' + key, 'm1'), value)

    @testing.requires_version('2.3.1')
    def test_undefined(self):
        cmd.fragment('gly')
        cmd.alter('name C', 'p["s_p_foo"] = "bar"')
        cmd.alter('name N', 'p["r_p_bar"] = 1.234')

        # save/load round-trip
        with testing.mktemp('.mae') as filename:
            cmd.save(filename)
            cmd.delete('*')
            cmd.load(filename, 'm1')

        myprops = {}
        cmd.iterate('all', 'myprops[name] = (p.s_p_foo, p.r_p_bar)', space=locals())

        self.assertEqual(myprops['O'], (None, None))
        self.assertEqual(myprops['C'], ("bar", None))
        self.assertEqual(myprops['N'], (None, 1.234))
