7
from pymol import cmd as global_cmd
import pymol


class Cmd:

    def __init__(self, _pymol, _COb):

        # store parent
        
        self._pymol = _pymol

        # store C object for easy access
    
        self._COb = _COb
        
        # private data
    
        self.color_sc = None
        self.reaper = None
        
        # CONSTANTS (pymol/constants.py)

        self.DEFAULT_ERROR = global_cmd.DEFAULT_ERROR
        self.DEFAULT_SUCCESS = global_cmd.DEFAULT_SUCCESS
        self.QuietException = global_cmd.QuietException
        self.Shortcut = global_cmd.Shortcut
        self._load2str = global_cmd._load2str
        self.boolean_dict = global_cmd.boolean_dict
        self.boolean_sc = global_cmd.boolean_sc
        self.fb_action = global_cmd.fb_action
        self.fb_mask = global_cmd.fb_mask
        self.fb_module = global_cmd.fb_module
        self.fb_debug = global_cmd.fb_debug # this cannot be right...
        self.file_ext_re = global_cmd.file_ext_re
        self.gz_ext_re = global_cmd.gz_ext_re
        self.loadable = global_cmd.loadable
        self.nt_hidden_path_re = global_cmd.nt_hidden_path_re
        self.palette_dict = global_cmd.palette_dict
        self.palette_sc = global_cmd.palette_sc
        self.parsing = global_cmd.parsing
        self.quote_alpha_list_re = global_cmd.quote_alpha_list_re
        self.repres = global_cmd.repres
        self.repres_sc = global_cmd.repres_sc
        self.safe_alpha_list_eval = global_cmd.safe_alpha_list_eval
        self.safe_list_eval = global_cmd.safe_list_eval
        self.safe_oname_re = global_cmd.safe_oname_re
        self.sanitize_alpha_list_re = global_cmd.sanitize_alpha_list_re
        self.sanitize_list_re = global_cmd.sanitize_list_re
        self.space_sc = global_cmd.space_sc
        self.stereo_dict = global_cmd.stereo_dict
        self.stereo_sc = global_cmd.stereo_sc
        self.toggle_dict = global_cmd.toggle_dict
        self.toggle_sc = global_cmd.toggle_sc
        self.window_dict = global_cmd.window_dict
        self.window_sc = global_cmd.window_sc

        
        # GLOBAL FUNCTIONS

        self.exp_path = global_cmd.exp_path

        # deferred initiailization
        
        global_cmd._deferred_init_pymol_internals(_pymol)
        
        # PRIVATE FUNCTIONS (requiring '_self' as a keyword argument)
        
        # locking.py

        self.reaper = None

        if 1:
            # use own locks (for performance)
            self.lock_api = _pymol.lock_api
            self.lock_api_c = _pymol.lock_api_c
            self.lock_api_data = _pymol.lock_api_data            
            self.lock_api_glut = _pymol.lock_api_glut
            self.lock_api_status = _pymol.lock_api_status
        else:
            # use global locks (for debugging)
            self.lock_api = global_cmd._pymol.lock_api
            self.lock_api_c = global_cmd._pymol.lock_api_c
            self.lock_api_data = global_cmd._pymol.lock_api_data            
            self.lock_api_glut = global_cmd._pymol.lock_api_glut
            self.lock_api_status = global_cmd._pymol.lock_api_status

        self.lock_api_allow_flush = 1

        self.lock_c = global_cmd.lock_c
        self.unlock_c = global_cmd.unlock_c
        self.lock_data = global_cmd.lock_data
        self.unlock_data = global_cmd.unlock_data
        self.lock_status_attempt = global_cmd.lock_status_attempt
        self.lock_status = global_cmd.lock_status
        self.unlock_status = global_cmd.unlock_status
        self.lock_glut = global_cmd.lock_glut
        self.unlock_glut = global_cmd.unlock_glut
        self.lock_without_glut = global_cmd.lock_without_glut
        self.lock = global_cmd.lock
        self.lock_attempt = global_cmd.lock_attempt
        self.unlock = global_cmd.unlock
        self.is_glut_thread = global_cmd.is_glut_thread
        self.setup_global_locks = global_cmd.setup_global_locks
        self.interrupt = global_cmd.interrupt
        self.block_flush = global_cmd.block_flush
        self.unblock_flush = global_cmd.unblock_flush        
        
        # pymol/checking.py

        self._raising = global_cmd._raising
        self.is_dict = global_cmd.is_dict
        self.is_error = global_cmd.is_error
        self.is_list = global_cmd.is_list
        self.is_ok = global_cmd.is_ok
        self.is_sequence = global_cmd.is_sequence
        self.is_string = global_cmd.is_string
        self.is_tuple = global_cmd.is_tuple

        # from pymol/feedingback.py

        self._feedback = global_cmd._feedback
        
        # from pymol/internal.py
        
        self._adjust_coord = global_cmd._adjust_coord
        self._coordset_update_spawn = global_cmd._coordset_update_spawn
        self._coordset_update_thread = global_cmd._coordset_update_thread
        self._copy_image = global_cmd._copy_image
        self._do = global_cmd._do
        self._dump_floats = global_cmd._dump_floats
        self._dump_ufloats = global_cmd._dump_ufloats
        self._fake_drag = global_cmd._fake_drag
        self._get_color_sc = global_cmd._get_color_sc
        self._get_feedback = global_cmd._get_feedback
        self._interpret_color = global_cmd._interpret_color
        self._interpret_color = global_cmd._interpret_color
        self._invalidate_color_sc = global_cmd._invalidate_color_sc
        self._invalidate_color_sc = global_cmd._invalidate_color_sc
        self._load = global_cmd._load
        self._mpng = global_cmd._mpng
        self._object_update_spawn = global_cmd._object_update_spawn
        self._object_update_thread = global_cmd._object_update_thread
        self._png = global_cmd._png
        self._quit = global_cmd._quit
        self._ray_anti_spawn = global_cmd._ray_anti_spawn
        self._ray_hash_spawn = global_cmd._ray_hash_spawn
        self._ray_spawn = global_cmd._ray_spawn
        self._refresh = global_cmd._refresh
        self._sgi_stereo = global_cmd._sgi_stereo
        self._validate_color_sc = global_cmd._validate_color_sc
        self._cache_set = global_cmd._cache_set
        self._cache_get = global_cmd._cache_get
        self._cache_clear = global_cmd._cache_clear
        self._cache_mark = global_cmd._cache_mark
        self._cache_purge = global_cmd._cache_purge

        # now we create the command langauge

        from pymol import keywords

        self.keyword = keywords.get_command_keywords()
        self.kw_list = self.keyword.keys()

        keywords.fix_list(self.kw_list)
        self.kwhash = self.Shortcut(self.kw_list)
        keywords.fix_dict(self.keyword)

        self.controlling = pymol.controlling
        self.completing = pymol.completing
        self.controlling = pymol.controlling
        self.editing = pymol.editing
        self.exporting = pymol.exporting
        self.moving = pymol.moving
        self.creating = pymol.creating
        self.viewing = pymol.viewing
        self.setting = pymol.setting
        self.commanding = pymol.commanding

        self.help_only = keywords.get_help_only_keywords()
        self.help_sc = self.Shortcut(self.keyword.keys()+self.help_only.keys())

        
        self.selection_sc = lambda sc=self.Shortcut,gn=self.get_names:sc(gn('public')+['all'])
        self.object_sc = lambda sc=self.Shortcut,gn=self.get_names:sc(gn('objects'))
        self.map_sc = lambda sc=self.Shortcut,gnot=self.get_names_of_type:sc(gnot('object:map'))
        self.contour_sc =  lambda sc=self.Shortcut,gnot=self.get_names_of_type:sc(
            gnot('object:mesh')+gnot('object:surface'))
        
        self.fb_action_sc = pymol.feedingback.fb_action_sc
        self.fb_module_sc = pymol.feedingback.fb_module_sc
        self.fb_mask_sc = pymol.feedingback.fb_mask_sc
        
        self.auto_arg = pymol.completing.get_auto_arg_list(self)
        self.color_sc = None

        # keyboard configuration
                
        from pymol import keyboard
        
        self.special = keyboard.get_special(self)

        self.shft_special = keyboard.get_shft_special(self)        
        self.alt_special = keyboard.get_alt_special(self)        
        self.ctrl_special = keyboard.get_ctrl_special(self)
        self.ctsh_special = keyboard.get_ctsh_special(self)

        self.ctrl = keyboard.get_ctrl(self)        
        self.alt = keyboard.get_alt(self)
        
# PUBLIC API METHODS which expect "self" as the first argument

# ========= WARNING WARNING WARNING WARNING ===========
# ========= AUTOGENERATED BEYOND THIS POINT ===========

    def _alt(self, *a, **k):
        k['_self']=self
        return apply(global_cmd._alt, a, k)
    
    def _ctrl(self, *a, **k):
        k['_self']=self
        return apply(global_cmd._ctrl, a, k)
    
    def _feedback(self, *a, **k):
        k['_self']=self
        return apply(global_cmd._feedback, a, k)
    
    def _special(self, *a, **k):
        k['_self']=self
        return apply(global_cmd._special, a, k)
    
    def abort(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.abort, a, k)
    
    def accept(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.accept, a, k)
    
    def alias(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.alias, a, k)
    
    def align(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.align, a, k)
    
    def alter(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.alter, a, k)
    
    def alter_list(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.alter_list, a, k)
    
    def alter_state(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.alter_state, a, k)
    
    def angle(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.angle, a, k)
    
    def attach(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.attach, a, k)
    
    def auto_measure(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.auto_measure, a, k)
    
    def backward(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.backward, a, k)
    
    def bg_color(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.bg_color, a, k)
    
    def bg_colour(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.bg_colour, a, k)
    
    def bond(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.bond, a, k)
    
    def button(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.button, a, k)
    
    def cache(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.cache, a, k)
    
    def capture(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.capture, a, k)
    
    def cartoon(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.cartoon, a, k)
    
    def cd(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.cd, a, k)
    
    def center(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.center, a, k)
    
    def check(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.check, a, k)
    
    def clean(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.clean, a, k)
    
    def clip(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.clip, a, k)
    
    def cls(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.cls, a, k)
    
    def color(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.color, a, k)
    
    def colour(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.colour, a, k)
    
    def commands(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.commands, a, k)
    
    def config_mouse(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.config_mouse, a, k)
    
    def copy(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.copy, a, k)
    
    def copy_image(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.copy_image, a, k)
    
    def count_atoms(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.count_atoms, a, k)
    
    def count_frames(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.count_frames, a, k)
    
    def count_states(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.count_states, a, k)
    
    def create(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.create, a, k)
    
    def cycle_valence(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.cycle_valence, a, k)
    
    def decline(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.decline, a, k)
    
    def del_colorection(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.del_colorection, a, k)
    
    def delete(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.delete, a, k)
    
    def deprotect(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.deprotect, a, k)
    
    def deselect(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.deselect, a, k)
    
    def dihedral(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dihedral, a, k)
    
    def dir(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dir, a, k)
    
    def dirty(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dirty, a, k)
    
    def dirty_wizard(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dirty_wizard, a, k)
    
    def disable(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.disable, a, k)
    
    def dist(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dist, a, k)
    
    def distance(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.distance, a, k)
    
    def do(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.do, a, k)
    
    def drag(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.drag, a, k)
    
    def draw(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.draw, a, k)
    
    def dss(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dss, a, k)
    
    def dummy(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dummy, a, k)
    
    def dump(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.dump, a, k)
    
    def edit(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.edit, a, k)
    
    def edit_mode(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.edit_mode, a, k)
    
    def enable(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.enable, a, k)
    
    def ending(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.ending, a, k)
    
    def export_coords(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.export_coords, a, k)
    
    def export_dots(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.export_dots, a, k)
    
    def extend(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.extend, a, k)
    
    def extract(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.extract, a, k)
    
    def fab(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.fab, a, k)
    
    def feedback(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.feedback, a, k)
    
    def fetch(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.fetch, a, k)
    
    def find_pairs(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.find_pairs, a, k)
    
    def finish_object(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.finish_object, a, k)
    
    def fit(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.fit, a, k)
    
    def fix_chemistry(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.fix_chemistry, a, k)
    
    def flag(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.flag, a, k)
    
    def forward(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.forward, a, k)
    
    def fragment(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.fragment, a, k)
    
    def frame(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.frame, a, k)
    
    def full_screen(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.full_screen, a, k)
    
    def fuse(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.fuse, a, k)
    
    def get(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get, a, k)
    
    def get_angle(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_angle, a, k)
    
    def get_area(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_area, a, k)
    
    def get_atom_coords(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_atom_coords, a, k)
    
    def get_bond_print(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_bond_print, a, k)
    
    def get_busy(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_busy, a, k)
    
    def get_chains(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_chains, a, k)
    
    def get_color_index(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_color_index, a, k)
    
    def get_color_indices(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_color_indices, a, k)
    
    def get_color_tuple(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_color_tuple, a, k)
    
    def get_colorection(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_colorection, a, k)
    
    def get_dihedral(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_dihedral, a, k)
    
    def get_distance(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_distance, a, k)
    
    def get_editor_scheme(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_editor_scheme, a, k)
    
    def get_extent(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_extent, a, k)
    
    def get_fastastr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_fastastr, a, k)
    
    def get_frame(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_frame, a, k)
    
    def get_idtf(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_idtf, a, k)
    
    def get_legal_name(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_legal_name, a, k)
    
    def get_modal_draw(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_modal_draw, a, k)
    
    def get_model(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_model, a, k)
    
    def get_movie_length(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_movie_length, a, k)
    
    def get_movie_locked(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_movie_locked, a, k)
    
    def get_movie_playing(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_movie_playing, a, k)
    
    def get_mtl_obj(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_mtl_obj, a, k)
    
    def get_names(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_names, a, k)
    
    def get_names_of_type(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_names_of_type, a, k)
    
    def get_object_color_index(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_object_color_index, a, k)
    
    def get_object_list(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_object_list, a, k)
    
    def get_object_matrix(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_object_matrix, a, k)
    
    def get_pdbstr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_pdbstr, a, k)
    
    def get_phipsi(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_phipsi, a, k)
    
    def get_position(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_position, a, k)
    
    def get_povray(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_povray, a, k)
    
    def get_progress(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_progress, a, k)
    
    def get_raw_alignment(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_raw_alignment, a, k)
    
    def get_renderer(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_renderer, a, k)
    
    def get_scene_dict(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_scene_dict, a, k)
    
    def get_scene_list(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_scene_list, a, k)
    
    def get_session(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_session, a, k)
    
    def get_setting_boolean(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_boolean, a, k)
    
    def get_setting_float(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_float, a, k)
    
    def get_setting_int(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_int, a, k)
    
    def get_setting_legacy(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_legacy, a, k)
    
    def get_setting_text(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_text, a, k)
    
    def get_setting_tuple(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_tuple, a, k)
    
    def get_setting_updates(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_setting_updates, a, k)
    
    def get_state(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_state, a, k)
    
    def get_symmetry(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_symmetry, a, k)
    
    def get_title(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_title, a, k)
    
    def get_type(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_type, a, k)
    
    def get_unused_name(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_unused_name, a, k)
    
    def get_version(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_version, a, k)
    
    def get_view(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_view, a, k)
    
    def get_vis(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_vis, a, k)
    
    def get_vrml(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_vrml, a, k)
    
    def get_wizard(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_wizard, a, k)
    
    def get_wizard_stack(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.get_wizard_stack, a, k)
    
    def gradient(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.gradient, a, k)
    
    def group(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.group, a, k)
    
    def h_add(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.h_add, a, k)
    
    def h_fill(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.h_fill, a, k)
    
    def h_fix(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.h_fix, a, k)
    
    def help(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.help, a, k)
    
    def hide(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.hide, a, k)
    
    def id_atom(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.id_atom, a, k)
    
    def identify(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.identify, a, k)
    
    def import_coords(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.import_coords, a, k)
    
    def index(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.index, a, k)
    
    def indicate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.indicate, a, k)
    
    def interrupt(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.interrupt, a, k)
    
    def intra_fit(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.intra_fit, a, k)
    
    def intra_rms(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.intra_rms, a, k)
    
    def intra_rms_cur(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.intra_rms_cur, a, k)
    
    def invert(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.invert, a, k)
    
    def isodot(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.isodot, a, k)
    
    def isolevel(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.isolevel, a, k)
    
    def isomesh(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.isomesh, a, k)
    
    def isosurface(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.isosurface, a, k)
    
    def iterate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.iterate, a, k)
    
    def iterate_state(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.iterate_state, a, k)
    
    def label(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.label, a, k)
    
    def label2(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.label2, a, k)
    
    def load(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load, a, k)
    
    def load_brick(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_brick, a, k)
    
    def load_callback(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_callback, a, k)
    
    def load_cgo(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_cgo, a, k)
    
    def load_coords(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_coords, a, k)
    
    def load_embedded(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_embedded, a, k)
    
    def load_map(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_map, a, k)
    
    def load_model(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_model, a, k)
    
    def load_object(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_object, a, k)
    
    def load_png(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_png, a, k)
    
    def load_raw(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_raw, a, k)
    
    def load_traj(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.load_traj, a, k)
    
    def log(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.log, a, k)
    
    def log_close(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.log_close, a, k)
    
    def log_open(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.log_open, a, k)
    
    def ls(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.ls, a, k)
    
    def madd(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.madd, a, k)
    
    def map_double(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.map_double, a, k)
    
    def map_halve(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.map_halve, a, k)
    
    def map_new(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.map_new, a, k)
    
    def map_set(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.map_set, a, k)
    
    def map_set_border(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.map_set_border, a, k)
    
    def map_trim(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.map_trim, a, k)
    
    def mappend(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mappend, a, k)
    
    def mask(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mask, a, k)
    
    def matrix_copy(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.matrix_copy, a, k)
    
    def matrix_reset(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.matrix_reset, a, k)
    
    def matrix_transfer(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.matrix_transfer, a, k)
    
    def mclear(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mclear, a, k)
    
    def mdo(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mdo, a, k)
    
    def mdump(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mdump, a, k)
    
    def mem(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mem, a, k)
    
    def meter_reset(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.meter_reset, a, k)
    
    def middle(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.middle, a, k)
    
    def mmatrix(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mmatrix, a, k)
    
    def mouse(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mouse, a, k)
    
    def move(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.move, a, k)
    
    def mplay(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mplay, a, k)
    
    def mpng(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mpng, a, k)
    
    def mray(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mray, a, k)
    
    def mset(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mset, a, k)
    
    def mstop(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mstop, a, k)
    
    def mtoggle(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mtoggle, a, k)
    
    def multisave(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.multisave, a, k)
    
    def mview(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.mview, a, k)
    
    def order(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.order, a, k)
    
    def orient(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.orient, a, k)
    
    def origin(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.origin, a, k)
    
    def overlap(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.overlap, a, k)
    
    def pair_fit(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.pair_fit, a, k)
    
    def phi_psi(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.phi_psi, a, k)
    
    def png(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.png, a, k)
    
    def pop(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.pop, a, k)
    
    def protect(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.protect, a, k)
    
    def pseudoatom(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.pseudoatom, a, k)
    
    def push_undo(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.push_undo, a, k)
    
    def pwd(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.pwd, a, k)
    
    def quit(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.quit, a, k)
    
    def ramp_new(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.ramp_new, a, k)
    
    def ray(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.ray, a, k)
    
    def read_mmodstr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.read_mmodstr, a, k)
    
    def read_molstr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.read_molstr, a, k)
    
    def read_pdbstr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.read_pdbstr, a, k)
    
    def read_sdfstr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.read_sdfstr, a, k)
    
    def read_xplorstr(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.read_xplorstr, a, k)
    
    def ready(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.ready, a, k)
    
    def rebuild(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rebuild, a, k)
    
    def recolor(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.recolor, a, k)
    
    def recolour(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.recolour, a, k)
    
    def redo(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.redo, a, k)
    
    def reference(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.reference, a, k)
    
    def refresh(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.refresh, a, k)
    
    def refresh_wizard(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.refresh_wizard, a, k)
    
    def reinitialize(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.reinitialize, a, k)
    
    def remove(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.remove, a, k)
    
    def remove_picked(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.remove_picked, a, k)
    
    def rename(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rename, a, k)
    
    def replace(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.replace, a, k)
    
    def replace_wizard(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.replace_wizard, a, k)
    
    def reset(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.reset, a, k)
    
    def resume(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.resume, a, k)
    
    def rewind(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rewind, a, k)
    
    def rms(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rms, a, k)
    
    def rms_cur(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rms_cur, a, k)
    
    def rock(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rock, a, k)
    
    def rotate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.rotate, a, k)
    
    def save(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.save, a, k)
    
    def scene(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.scene, a, k)
    
    def scene_order(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.scene_order, a, k)
    
    def sculpt_activate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.sculpt_activate, a, k)
    
    def sculpt_deactivate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.sculpt_deactivate, a, k)
    
    def sculpt_iterate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.sculpt_iterate, a, k)
    
    def sculpt_purge(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.sculpt_purge, a, k)
    
    def select(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.select, a, k)
    
    def select_list(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.select_list, a, k)
    
    def set(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set, a, k)
    
    def set_bond(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_bond, a, k)
    
    def set_color(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_color, a, k)
    
    def set_colorection(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_colorection, a, k)
    
    def set_colorection_name(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_colorection_name, a, k)
    
    def set_colour(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_colour, a, k)
    
    def set_colour(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_colour, a, k)
    
    def set_dihedral(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_dihedral, a, k)
    
    def set_geometry(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_geometry, a, k)
    
    def set_key(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_key, a, k)
    
    def set_name(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_name, a, k)
    
    def set_object_color(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_object_color, a, k)
    
    def set_object_color(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_object_color, a, k)
    
    def set_object_ttt(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_object_ttt, a, k)
    
    def set_session(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_session, a, k)
    
    def set_symmetry(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_symmetry, a, k)
    
    def set_title(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_title, a, k)
    
    def set_view(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_view, a, k)
    
    def set_vis(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_vis, a, k)
    
    def set_wizard(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_wizard, a, k)
    
    def set_wizard_stack(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.set_wizard_stack, a, k)
    
    def show(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.show, a, k)
    
    def show_as(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.show_as, a, k)
    
    def show_help(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.show_help, a, k)
    
    def slice_new(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.slice_new, a, k)
    
    def smooth(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.smooth, a, k)
    
    def sort(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.sort, a, k)
    
    def space(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.space, a, k)
    
    def spectrum(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.spectrum, a, k)
    
    def spheroid(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.spheroid, a, k)
    
    def splash(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.splash, a, k)
    
    def split_states(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.split_states, a, k)
    
    def stereo(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.stereo, a, k)
    
    def super(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.super, a, k)
    
    def symexp(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.symexp, a, k)
    
    def sync(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.sync, a, k)
    
    def system(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.system, a, k)
    
    def test(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.test, a, k)
    
    def toggle(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.toggle, a, k)
    
    def torsion(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.torsion, a, k)
    
    def transform_object(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.transform_object, a, k)
    
    def transform_selection(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.transform_selection, a, k)
    
    def translate(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.translate, a, k)
    
    def translate_atom(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.translate_atom, a, k)
    
    def turn(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.turn, a, k)
    
    def unbond(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.unbond, a, k)
    
    def undo(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.undo, a, k)
    
    def ungroup(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.ungroup, a, k)
    
    def unmask(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.unmask, a, k)
    
    def unpick(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.unpick, a, k)
    
    def unprotect(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.unprotect, a, k)
    
    def unset(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.unset, a, k)
    
    def unset_bond(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.unset_bond, a, k)
    
    def update(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.update, a, k)
    
    def valence(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.valence, a, k)
    
    def vdw_fit(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.vdw_fit, a, k)
    
    def view(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.view, a, k)
    
    def viewport(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.viewport, a, k)
    
    def window(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.window, a, k)
    
    def wizard(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.wizard, a, k)
    
    def zoom(self, *a, **k):
        k['_self']=self
        return apply(global_cmd.zoom, a, k)
    
