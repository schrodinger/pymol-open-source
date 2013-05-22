'''
Testing: pymol.creating
'''

from pymol import cmd, testing, stored

import unittest

class TestCreating(testing.PyMOLTestCase):

    def testGroup(self):
        cmd.pseudoatom('m1')
        cmd.pseudoatom('m2')
        cmd.pseudoatom('m3')
        m_list = []
        cmd.group('g1', 'm1 m2')
        cmd.iterate('g1', 'm_list.append(model)', space=locals())
        self.assertItemsEqual(m_list, ['m1', 'm2'])
        m_list = []
        cmd.ungroup('m2')
        cmd.iterate('g1', 'm_list.append(model)', space=locals())
        self.assertItemsEqual(m_list, ['m1'])

    def testUngroup(self):
        # see testGroup
        pass

    @unittest.skip("Not yet implemented.")
    def testIsomesh(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testVolume(self):
        pass
    
    @unittest.skip("Not yet implemented.")
    def testIsosurface(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testIsodot(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testIsolevel(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testGradient(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testCopy(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testSymexp(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testFragment(self):
        pass

    def testCreate(self):
        cmd.fragment("ala")
        cmd.create("foo", "ala", 1, 1)

        names = cmd.get_names()
        names.sort()

        self.assertEquals(names, ["ala", "foo"])
        
        # TODO: Create a function to compare to chempy models
        # m1 = cmd.get_model("foo")
        # m2 = cmd.get_model("ala")
        # self.assertChempyModelEquals(m1, m2)

        self.assertEquals(cmd.fit("ala", "foo"), 0.0)

        self.assertEquals(cmd.count_states(), 1)

    def testCreateMultipleStates(self):
        cmd.fragment("ala")

        # simple creation

        self.assertEquals(cmd.count_states(), 1)

        cmd.create("foo", "ala", 1, 1)

        self.assertEquals(cmd.count_states(), 1)

        # creation into state 1

        cmd.create("foo", "ala", 1, 1)

        self.assertEquals(cmd.count_states(), 1)

        # creation into state 2

        cmd.create("foo", "ala", 1, 2)

        self.assertEquals(cmd.count_states(), 2)
        
    def create_many_states(self, fragment, obj, count):
        cmd.fragment(fragment)
        for x in range(count):
            cmd.create(obj, fragment, 1, count)

    def testCreateMany(self):

        # 10 states
        self.create_many_states("thr", "foo", 10)

        self.assertEquals(cmd.count_states("foo"), 10)
        self.assertEquals(cmd.count_states(), 10)

        # 100 states
        self.create_many_states("trp", "bar", 100)

        self.assertEquals(cmd.count_states("foo"), 10)
        self.assertEquals(cmd.count_states("bar"), 100)
        self.assertEquals(cmd.count_states(), 100)

        # 1000 states
        
        self.create_many_states("ile", "bam", 1000)

        self.assertEquals(cmd.count_states(), 1000)
        

    def testCreateLarge(self):
        pass

    def testCreateManyLarge(self):
        pass
        


    @unittest.skip("Not yet implemented.")
    def testExtract(self):
        pass
    
    @unittest.skip("Not yet implemented.")
    def testUnquote(self):
        pass

    @unittest.skip("Not yet implemented.")
    def testPseudoatom(self):
        pass

