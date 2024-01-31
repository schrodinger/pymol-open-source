'''
Testing setting/getting properties for different types
'''

import random
import unittest
from pymol import cmd, testing

max_threads = cmd.get_setting_int('max_threads')

@testing.requires('gui', 'no_run_all', 'multicore')
class TestAsyncBuilds(testing.PyMOLTestCase):

    @testing.foreach.product(['surface', 'cartoon'], [0, 1])
    def testAsyncBuilds(self, rep, async_builds):
        target = "1aon"
        if async_builds:
            msg = '%s cpus' % max_threads
        else:
            msg = '1 cpu'
        cmd.set("async_builds", async_builds)
    
        cmd.load(self.datafile('1aon.pdb.gz'), target)
        for x in cmd.get_chains():
            cmd.create("Chain_%s" % x, target + " & c. " + x)
        cmd.delete(target)

        with self.timing('%s' % msg):
            cmd.show_as(rep)
            cmd.draw()
