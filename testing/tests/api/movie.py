from pymol import movie, cmd, testing, stored

class TestMovie(testing.PyMOLTestCase):

    def _make_n_states(self, nstates):
        cmd.fragment('gly', 'm1')
        for state in range(2, nstates + 1):
            cmd.create('m1', 'm1', 1, state)

    def test_get_movie_fps(self):
        self.assertEqual(30, movie.get_movie_fps(_self=cmd))
        cmd.set('movie_fps', 10) # > 0
        self.assertEqual(10, movie.get_movie_fps(_self=cmd))
        cmd.set('movie_fps', 0) # <= 0
        self.assertEqual(30, movie.get_movie_fps(_self=cmd))

    def test_sweep(self):
        N = 3
        self._make_n_states(N)
        movie.sweep()
        self.assertEqual(cmd.count_frames(), N * 2)
        P = 5
        movie.sweep(pause=P)
        self.assertEqual(cmd.count_frames(), N * 2 + P * 2)
        C = 3
        movie.sweep(cycles=C)
        self.assertEqual(cmd.count_frames(), N * 2 * C)

    def test_pause(self):
        N = 3
        self._make_n_states(N)
        P = 10
        movie.pause(pause=P)
        self.assertEqual(cmd.count_frames(), 2 * P + N)

    def test_rock(self):
        self.skipTest("TODO") # movie.rock(first=1,last=-1,angle=30,phase=0,loop=1,axis='y'):

    def test_roll(self):
        self.skipTest("TODO") # movie.roll(first=1,last=-1,loop=1,axis='y'):

    def test_tdroll(self):
        self.skipTest("TODO") # movie.tdroll(first,rangex,rangey,rangez,skip=1):

    def test_zoom(self):
        self.skipTest("TODO") # movie.zoom(first,last,step=1,loop=1,axis='z'):

    def test_nutate(self):
        self.skipTest("TODO") # movie.nutate(first,last,angle=30,phase=0,loop=1,shift=math.pi/2.0):

    def test_screw(self):
        self.skipTest("TODO") # movie.screw(first,last,step=1,angle=30,phase=0,loop=1,axis='y'):

    def test_timed_roll(self):
        self.skipTest("TODO") # movie.timed_roll(period=12.0,cycles=1,axis='y'):

    def test_add_blank(self):
        self.skipTest("TODO") # movie.add_blank(duration=12.0,start=0):

    def test_add_roll(self):
        self.skipTest("TODO") # movie.add_roll(duration=12.0,loop=1,axis='y',start=0):

    def test_add_rock(self):
        self.skipTest("TODO") # movie.add_rock(duration=8.0,angle=30.0,loop=1,axis='y',start=0):

    def test_add_state_sweep(self):
        self.skipTest("TODO") # movie.add_state_sweep(factor=1,pause=2.0,first=-1,last=-1,loop=1,start=0):

    def test_add_state_loop(self):
        self.skipTest("TODO") # movie.add_state_loop(factor=1,pause=2.0,first=-1,last=-1,loop=1,start=0):

    def test_add_nutate(self):
        self.skipTest("TODO") # movie.add_nutate(duration=8.0, angle=30.0, spiral=0, loop=1, 

    @testing.requires_version('2.4')
    def test_add_scenes(self):
        cmd.fragment('gly', 'm1')
        cmd.scene('001', 'store')
        cmd.turn('x', 90)
        cmd.scene('002', 'store')
        movie.add_scenes()
        self.assertEqual(cmd.count_frames(), 615)
        cmd.mset()
        movie.add_scenes(['001', '002'], pause=3)
        self.assertEqual(cmd.count_frames(), 315)
        cmd.mset()
        movie.add_scenes('["001", "002"]') # string
        self.assertEqual(cmd.count_frames(), 615)

    def test_produce(self):
        self.skipTest("TODO") # movie.produce(filename, mode='', first=0, last=0, preserve=0,
