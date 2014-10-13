'''
Test that all reps are present.

PYMOL-2222
'''

import os
from pymol import cmd, testing, stored

REPS = [
    'sticks',
    'labels',
    'surface',
    'ribbon',
    # 'cartoon',
    # 'ellipsoids',     # PYMOL-1356
    'dots',
    'nonbonded',
    'nb_spheres',
    'mesh',
    'lines',
    'spheres',

    'dashes',
    'angles',
    'dihedrals',

    # 'cell',
    # 'extent',

    # 'slice',
    # 'callback',
    # 'cgo',
    # 'volume',
]

CARTOONS = [
    'putty',
    'rectangle',
    'arrow',
    'oval',
    # 'skip',
    'tube',
    'dumbbell',
    # 'automatic',
    'loop',
]

class TestReps(testing.PyMOLTestCase):

    def testRepsExist(self):
        cmd.viewport(200, 150)
        cmd.load(self.datafile('1oky-frag.pdb'), 'm1')

        # make some nonbonded
        cmd.unbond('resi 115-', 'resi 115-')

        # labels
        cmd.label('all', 'name')

        # measurements
        cmd.distance('measure1', 'index 1', 'index 10')
        cmd.angle('measure1', 'index 1', 'index 10', 'index 20')
        cmd.dihedral('measure1', 'index 1', 'index 10', 'index 20', 'index 30')

        # color test setup
        cmd.color('white', '*')
        cmd.set('ambient', 1)
        cmd.set('depth_cue', 0)
        cmd.set('antialias', 0)
        cmd.set('line_smooth', 0)
        cmd.orient()

        # test most reps
        for rep in REPS:
            cmd.show_as(rep)
            self.assertImageHasColor('white', msg='rep missing: ' + rep)

        # test cartoon
        cmd.show_as('cartoon')
        for cart in CARTOONS:
            cmd.cartoon(cart)
            self.assertImageHasColor('white', msg='cartoon missing: ' + cart)
