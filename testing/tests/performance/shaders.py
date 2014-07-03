import os
import pymol
from pymol import cmd, testing, stored

@testing.requires("no_run_all", "gui")
class StressShaders(testing.PyMOLTestCase):

    def _setup_movie(self, rep):
        cmd.load(self.datafile('1aon.pdb.gz'))
        cmd.show_as(rep)
        cmd.orient()
        cmd.mset('1x100')
        pymol.movie.roll(1, cmd.count_frames(), 0, 'y')

    @testing.foreach(
            (0, 0),
            (0, 0, "surface"),
            (1, 0),
            (1, 0, "surface"),
            (1, 1),
            )
    def test(self, use_shaders, cgo_shader_ub, rep="cartoon"):
        self._setup_movie(rep)

        cmd.set("use_shaders", use_shaders)
        cmd.set("cgo_shader_ub_color", cgo_shader_ub)
        cmd.set("cgo_shader_ub_flags", cgo_shader_ub)
        cmd.set("cgo_shader_ub_normal", cgo_shader_ub)

        with testing.mkdtemp() as tempdir, self.timing():
            cmd.mpng(os.path.join(tempdir, "f"))
