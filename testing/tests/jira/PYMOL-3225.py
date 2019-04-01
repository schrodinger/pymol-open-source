'''
32bit saved PSE with negative atom flags (unsigned casted to signed)
'''

from pymol import cmd, testing

class Test3225(testing.PyMOLTestCase):

    @testing.foreach((32,), (64,))
    def test(self, bits):
        cmd.load(self.datafile('AC-1744-{}bit.pse.gz'.format(bits)))
        flags_list = []
        cmd.iterate(
            'name CA',
            'flags_list.append(flags & 0xFFFFFFFF)',
            space=locals())
        self.assertEqual(flags_list, [0x88000040] * 2)
