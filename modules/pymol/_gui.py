'''
PyMOL GUI Data (toolkit independant)
'''

import sys
import os
import webbrowser

class PyMOLDesktopGUI(object):
    '''Superclass for PyMOL Desktop Applications'''

    file_open = None
    file_autoload_mtz = None
    file_fetch_pdb = None
    session_save = None
    session_save_as = None
    file_save = None
    file_save_map = None
    file_save_aln = None
    file_save_png = None
    file_save_wrl = None
    file_save_dae = None
    file_save_pov = None
    file_save_mpeg = None
    file_save_mov = None
    file_save_mpng = None
    file_save_gltf = None
    file_save_stl = None
    log_open = None
    log_resume = None
    log_append = None
    file_run = None
    edit_pymolrc = None
    confirm_quit = None
    settings_edit_all_dialog = None
    edit_colors_dialog = None
    cd_dialog = None
    show_about = None

    def new_window(self, extra_argv=()):
        import pymol

        python = sys.executable
        if os.path.isfile(python + 'w'):
            # fixes menu focus on macOS
            python += 'w'

        args = [python, pymol.__file__, '-N', pymol.invocation.options.gui
                ] + list(extra_argv)

        os.spawnv(os.P_NOWAITO, args[0], args)

    def get_menudata(self, cmd=None):
        '''Get the top level application menu as a list data structure'''

        if cmd is None:
            cmd = self.cmd

        def F_scene_menu(action):
            return [
                ('command', k, lambda k=k, a=action: cmd.scene(k, a))
                for k in ['F' + str(i) for i in range(1, 13)]
            ]

        def transparency_menu(setting_name):
            return [
                ('radio', lab, setting_name, val)
                for lab, val in [ ('Off', 0.0), ('20%', 0.2), ('40%', 0.4),
                    ('50%', 0.5), ('60%', 0.6), ('80%', 0.8) ]
            ]

        file_browser_cmd = (
                'explorer .'    if sys.platform == 'win32' else
                'open .'        if sys.platform == 'darwin' else
                'xdg-open .')

        return [
            ('menu', 'File', [
                ('menu', 'New PyMOL Window',  [
                    ('command', 'Default',                  self.new_window),
                    ('command', 'Ignore .pymolrc and plugins (-k)', lambda: self.new_window(('-k',))),
                ]),
                ('separator',),
                ('command', 'Open...',                      self.file_open),
                ('open_recent_menu',),
                ('command', 'Get PDB...',                   self.file_fetch_pdb),
                ('separator',),
                ('command', 'Save Session',                 self.session_save),
                ('command', 'Save Session As...',           self.session_save_as),
                ('separator',),
                ('command', 'Export Molecule...',             self.file_save),
                ('command', 'Export Map...',                  self.file_save_map),
                ('command', 'Export Alignment...',            self.file_save_aln),
                ('menu', 'Export Image As', [
                    ('command', 'PNG...',           self.file_save_png),
                    ('separator',),
                    ('command', 'VRML 2...',        self.file_save_wrl),
                    ('command', 'COLLADA...',       self.file_save_dae),
                    ('command', 'GLTF...',          self.file_save_gltf),
                    ('command', 'POV-Ray...',       self.file_save_pov),
                    ('command', 'STL...',           self.file_save_stl),
                ]),
                ('menu', 'Export Movie As', [
                    ('command', 'MPEG...',          self.file_save_mpeg),
                    ('command', 'Quicktime...',     self.file_save_mov),
                    ('separator',),
                    ('command', 'PNG Images...',    self.file_save_mpng),
                ]),
                ('separator',),
                ('menu', 'Log File', [
                    ('command', 'Open...',          self.log_open),
                    ('command', 'Resume...',        self.log_resume),
                    ('command', 'Append...',        self.log_append),
                    ('command', 'Close',            cmd.log_close),
                ]),
                ('command', 'Run Script...',        self.file_run),
                ('menu', 'Working Directory', [
                    ('command', 'Change...',        self.cd_dialog),
                    ('command', 'File Browser',     lambda: self.cmd.system(file_browser_cmd)),
                ]),
                ('separator',),
                ('command', 'Edit pymolrc',         self.edit_pymolrc),
                ('separator',),
                ('menu', 'Reinitialize', [
                    ('command', 'Everything',           cmd.reinitialize),
                    ('command', 'Original Settings',    'reinitialize original_settings'),
                    ('command', 'Stored Settings',      'reinitialize settings'),
                    ('separator',),
                    ('command', 'Store Current Settings',   'reinitialize store_defaults'),
                ]),
                ('command', 'Quit',                 self.confirm_quit),
            ]),
            ('menu', 'Edit', [
                ('command', 'Undo [Ctrl-Z]', cmd.undo),
                ('command', 'Redo [Ctrl-Y]', cmd.redo),
            ]),
            ('menu', 'Build', [
                ('menu', 'Fragment', [
                    ('command', 'Acetylene [Alt-J]', "editor.attach_fragment('pk1','acetylene',2,0)"),
                    ('command', 'Amide N->C [Alt-1]', "editor.attach_fragment('pk1','formamide',3,1)"),
                    ('command', 'Amide C->N [Alt-2]', "editor.attach_fragment('pk1','formamide',5,0)"),
                    ('command', 'Bromine [Ctrl-Shift-B]', "replace Br,1,1"),
                    ('command', 'Carbon [Ctrl-Shift-C]', "replace C,4,4"),
                    ('command', 'Carbonyl [Alt-0]', "editor.attach_fragment('pk1','formaldehyde',2,0)"),
                    ('command', 'Chlorine [Ctrl-Shift-L]', "replace Cl,1,1"),
                    ('command', 'Cyclobutyl [Alt-4]', "editor.attach_fragment('pk1','cyclobutane',4,0)"),
                    ('command', 'Cyclopentyl [Alt-5]', "editor.attach_fragment('pk1','cyclopentane',5,0)"),
                    ('command', 'Cyclopentadiene [Alt-8]', "editor.attach_fragment('pk1','cyclopentadiene',5,0)"),
                    ('command', 'Cyclohexyl [Alt-6]', "editor.attach_fragment('pk1','cyclohexane',7,0)"),
                    ('command', 'Cycloheptyl [Alt-7]', "editor.attach_fragment('pk1','cycloheptane',8,0)"),
                    ('command', 'Fluorine [Ctrl-Shift-F]', "replace F,1,1"),
                    ('command', 'Iodine [Ctrl-Shift-I]', "replace I,1,1"),
                    ('command', 'Methane [Ctrl-Shift-M]', "editor.attach_fragment('pk1','methane',1,0)"),
                    ('command', 'Nitrogen [Ctrl-Shift-N]', "replace N,4,3"),
                    ('command', 'Oxygen [Ctrl-Shift-O]', "replace O,4,2"),
                    ('command', 'Sulfer [Ctrl-Shift-S]', "replace S,2,2"),
                    ('command', 'Sulfonyl [Alt-3]', "editor.attach_fragment('pk1','sulfone',3,1)"),
                    ('command', 'Phosphorus [Ctrl-Shift-P]', "replace P,4,3"),
                ]),
                ('menu', 'Residue', [
                    ('command', lab, lambda v=val: cmd.editor.attach_amino_acid('pk1', v))
                    for lab, val in [
                        ('Acetyl [Alt-B]', 'ace'),
                        ('Alanine [Alt-A]', 'ala'),
                        ('Amine', 'nhh'),
                        ('Aspartate [Alt-D]', 'asp'),
                        ('Asparagine [Alt-N]', 'asn'),
                        ('Arginine [Alt-R]', 'arg'),
                        ('Cysteine [Alt-C]', 'cys'),
                        ('Glutamate [Alt-E]', 'glu'),
                        ('Glutamine [Alt-Q]', 'gln'),
                        ('Glycine [Alt-G]', 'gly'),
                        ('Histidine [Alt-H]', 'his'),
                        ('Isoleucine [Alt-I]', 'ile'),
                        ('Leucine [Alt-L]', 'leu'),
                        ('Lysine [Alt-K]', 'lys'),
                        ('Methionine [Alt-M]', 'met'),
                        ('N-Methyl [Alt-Z]', 'nme'),
                        ('Phenylalanine [Alt-F]', 'phe'),
                        ('Proline [Alt-P]', 'pro'),
                        ('Serine [Alt-S]', 'ser'),
                        ('Threonine [Alt-T]', 'thr'),
                        ('Tryptophan [Alt-W]', 'trp'),
                        ('Tyrosine [Alt-Y]', 'tyr'),
                        ('Valine [Alt-V]', 'val'),
                    ]
                ] + [
                    ('separator',),
                ] + [
                    ('radio', 'Helix', 'secondary_structure', 1),
                    ('radio', 'Antiparallel Beta Sheet', 'secondary_structure', 2),
                    ('radio', 'Parallel Beta Sheet', 'secondary_structure', 3),
                ]),
                ('separator',),
                ('menu', 'Sculpting', [
                    ('check', 'Auto-Sculpting', 'auto_sculpt'),
                    ('check', 'Sculpting', 'sculpting'),
                    ('separator',),
                    ('command', 'Activate', 'sculpt_activate all'),
                    ('command', 'Deactivate', 'sculpt_deactivate all'),
                    ('command', 'Clear Memory', cmd.sculpt_purge),
                    ('separator',),
                    ('radio', '1 Cycle per Update', 'sculpting_cycles', 1),
                ] + [
                    ('radio', str(val) + ' Cycles per Update', 'sculpting_cycles', val)
                    for val in [3, 10, 33, 100, 333, 1000]
                ] + [
                    ('separator',),
                ] + [
                    ('radio', lab, 'sculpt_field_mask', val)
                    for lab, val in [
                        ('Bonds Only',              0x01),
                        ('Bonds and Angles Only',   0x01|0x02),
                        ('Local Geometry Only',     0x20 - 1),
                        ('All Except VDW',                  ~(0x20 | 0x40)),
                        ('All Except 1-4 VDW and Torsions', ~(0x40 | 0x80)),
                        ('All Terms',               0xFF),
                    ]
                ]),
                ('separator',),
                ('command', 'Cycle Bond Valence [Ctrl-Shift-W]', "cycle_valence"),
                ('command', 'Fill Hydrogens on (pk1) [Ctrl-Shift-R]', "h_fill"),
                ('command', 'Invert (pk2)-(pk1)-(pk3) [Ctrl-Shift-E]', "invert"),
                ('command', 'Create Bond (pk1)-(pk2) [Ctrl-Shift-T]', "bond"),
                ('separator',),
                ('command', 'Remove (pk1) [Ctrl-Shift-D]', 'remove pk1'),
                ('separator',),
                ('command', 'Make (pk1) Positive [Ctrl-Shift-K]',   'alter pk1, formal_charge=1.'),
                ('command', 'Make (pk1) Negative [Ctrl-Shift-J]',   'alter pk1, formal_charge=-1.'),
                ('command', 'Make (pk1) Neutral [Ctrl-Shift-U]',    'alter pk1, formal_charge=0.'),
            ]),
            ('menu', 'Movie', [
                ('menu', 'Append', [
                    ('command', str(i) + ' second', lambda i=i: cmd.movie.add_blank(i))
                    for i in [.25, .5, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 48, 60]
                ]),
                ('separator',),
                ('menu', 'Program', [
                    ('menu', 'Camera Loop', [
                        ('menu', 'Nutate', [
                            ('command', '15 deg. over 4 sec.', lambda: self.mvprg("movie.add_nutate(4,15,start=%d)")),
                            ('command', '15 deg. over 8 sec.', lambda: self.mvprg("movie.add_nutate(8,15,start=%d)")),
                            ('command', '15 deg. over 12 sec.', lambda: self.mvprg("movie.add_nutate(12,15,start=%d)")),
                            ('separator', ''),
                            ('command', '30 deg. over 4 sec.', lambda: self.mvprg("movie.add_nutate(4,30,start=%d)")),
                            ('command', '30 deg. over 8 sec.', lambda: self.mvprg("movie.add_nutate(8,30,start=%d)")),
                            ('command', '30 deg. over 12 sec.', lambda: self.mvprg("movie.add_nutate(12,30,start=%d)")),
                            ('command', '30 deg. over 16 sec.', lambda: self.mvprg("movie.add_nutate(16,30,start=%d)")),
                            ('separator', ''),
                            ('command', '60 deg. over 8 sec.', lambda: self.mvprg("movie.add_nutate(8,60,start=%d)")),
                            ('command', '60 deg. over 16 sec.', lambda: self.mvprg("movie.add_nutate(16,60,start=%d)")),
                            ('command', '60 deg. over 24 sec.', lambda: self.mvprg("movie.add_nutate(24,60,start=%d)")),
                            ('command', '60 deg. over 32 sec.', lambda: self.mvprg("movie.add_nutate(32,60,start=%d)")),
                        ]),
                        ('separator',),
                        ('menu', 'X-Rock', [
                            ('command', '30 deg. over 2 sec.', lambda: self.mvprg("movie.add_rock(2,30,axis='x',start=%d)")),
                            ('command', '30 deg. over 4 sec.', lambda: self.mvprg("movie.add_rock(4,30,axis='x',start=%d)")),
                            ('command', '30 deg. over 8 sec.', lambda: self.mvprg("movie.add_rock(8,30,axis='x',start=%d)")),
                            ('separator',),
                            ('command', '60 deg. over 4 sec.', lambda: self.mvprg("movie.add_rock(4,60,axis='x',start=%d)")),
                            ('command', '60 deg. over 8 sec.', lambda: self.mvprg("movie.add_rock(8,60,axis='x',start=%d)")),
                            ('command', '60 deg. over 16 sec.', lambda: self.mvprg("movie.add_rock(16,60,axis='x',start=%d)")),
                            ('separator',),
                            ('command', '90 deg. over 6 sec.', lambda: self.mvprg("movie.add_rock(6,90,axis='x',start=%d)")),
                            ('command', '90 deg. over 12 sec.', lambda: self.mvprg("movie.add_rock(12,90,axis='x',start=%d)")),
                            ('command', '90 deg. over 24 sec.', lambda: self.mvprg("movie.add_rock(24,90,axis='x',start=%d)")),
                            ('separator',),
                            ('command', '120 deg. over 8 sec.', lambda: self.mvprg("movie.add_rock(8,120,axis='x',start=%d)")),
                            ('command', '120 deg. over 16 sec.', lambda: self.mvprg("movie.add_rock(16,120,axis='x',start=%d)")),
                            ('command', '120 deg. over 32 sec.', lambda: self.mvprg("movie.add_rock(32,120,axis='x',start=%d)")),
                            ('separator',),
                            ('command', '180 deg. over 12 sec.', lambda: self.mvprg("movie.add_rock(12,179.99,axis='x',start=%d)")),
                            ('command', '180 deg. over 24 sec.', lambda: self.mvprg("movie.add_rock(24,179.99,axis='x',start=%d)")),
                            ('command', '180 deg. over 48 sec.', lambda: self.mvprg("movie.add_rock(48,179.99,axis='x',start=%d)")),
                        ]),
                        ('menu', 'X-Roll', [
                            ('command', '4 seconds', lambda: self.mvprg("movie.add_roll(4.0,axis='x',start=%d)")),
                            ('command', '8 seconds', lambda: self.mvprg("movie.add_roll(8.0,axis='x',start=%d)")),
                            ('command', '16 seconds', lambda: self.mvprg("movie.add_roll(16.0,axis='x',start=%d)")),
                            ('command', '32 seconds', lambda: self.mvprg("movie.add_roll(32.0,axis='x',start=%d)")),
                        ]),
                        ('separator',),
                        ('menu', 'Y-Rock', [
                            ('command', '30 deg. over 2 sec.', lambda: self.mvprg("movie.add_rock(2,30,axis='y',start=%d)")),
                            ('command', '30 deg. over 4 sec.', lambda: self.mvprg("movie.add_rock(4,30,axis='y',start=%d)")),
                            ('command', '30 deg. over 8 sec.', lambda: self.mvprg("movie.add_rock(8,30,axis='y',start=%d)")),
                            ('separator',),
                            ('command', '60 deg. over 4 sec.', lambda: self.mvprg("movie.add_rock(4,60,axis='y',start=%d)")),
                            ('command', '60 deg. over 8 sec.', lambda: self.mvprg("movie.add_rock(8,60,axis='y',start=%d)")),
                            ('command', '60 deg. over 16 sec.', lambda: self.mvprg("movie.add_rock(16,60,axis='y',start=%d)")),
                            ('separator',),
                            ('command', '90 deg. over 6 sec.', lambda: self.mvprg("movie.add_rock(6,90,axis='y',start=%d)")),
                            ('command', '90 deg. over 12 sec.', lambda: self.mvprg("movie.add_rock(12,90,axis='y',start=%d)")),
                            ('command', '90 deg. over 24 sec.', lambda: self.mvprg("movie.add_rock(24,90,axis='y',start=%d)")),
                            ('separator',),
                            ('command', '120 deg. over 8 sec.', lambda: self.mvprg("movie.add_rock(8,120,axis='y',start=%d)")),
                            ('command', '120 deg. over 16 sec.', lambda: self.mvprg("movie.add_rock(16,120,axis='y',start=%d)")),
                            ('command', '120 deg. over 32 sec.', lambda: self.mvprg("movie.add_rock(32,120,axis='y',start=%d)")),
                            ('separator',),
                            ('command', '180 deg. over 12 sec.', lambda: self.mvprg("movie.add_rock(12,179.99,axis='y',start=%d)")),
                            ('command', '180 deg. over 24 sec.', lambda: self.mvprg("movie.add_rock(24,179.99,axis='y',start=%d)")),
                            ('command', '180 deg. over 48 sec.', lambda: self.mvprg("movie.add_rock(48,179.99,axis='y',start=%d)")),
                            ]),
                        ('menu', 'Y-Roll', [
                            ('command', '4 seconds', lambda: self.mvprg("movie.add_roll(4.0,axis='y',start=%d)")),
                            ('command', '8 seconds', lambda: self.mvprg("movie.add_roll(8.0,axis='y',start=%d)")),
                            ('command', '16 seconds', lambda: self.mvprg("movie.add_roll(16.0,axis='y',start=%d)")),
                            ('command', '32 seconds', lambda: self.mvprg("movie.add_roll(32.0,axis='y',start=%d)")),
                        ]),
                    ]),
                    ('separator',),
                    ('menu', 'Scene Loop', [
                        ('menu', lab, [
                            ('command', '%d deg. over %d sec.' % (angle, sec), lambda
                                s='set sweep_angle,%d;cmd.movie.add_scenes(None, %d, rock=%d, start=%%d)' % (angle, sec, rock):
                                self.mvprg(s))
                            for angle, seconds in ((30, (2,4,8)), (60, (4,8,16)), (90, (6,12,24)), (120, (8,16,32)))
                            for sec in seconds
                        ])
                        for lab, rock in [
                            ('Nutate', 4),
                            ('X-Rock', 2),
                            ('Y-Rock', 1),
                        ]
                    ] + [
                        ('menu', 'Steady', [
                            ('command',  str(val) + ' seconds each',
                                lambda s="movie.add_scenes(None,%.1f," % val: self.mvprg(s + "rock=0,start=%d)"))
                            for val in [1, 2, 4, 8, 12, 16, 24]
                        ]),
                    ]),
                    ('separator',),
                ] + [
                    ('menu', lab, [
                        ('menu', "Full Speed" if speed == 1 else ("1/%d Speed" % speed), [
                            ('command', (str(pause) + " second pause") if pause else "no pause",
                                lambda s=fmt % (speed, pause): self.mvprg(s + ", start=%d)"))
                            for pause in [0, 1, 2, 4]
                            ])
                        for speed in [1, 2, 3, 4, 8, 16]
                    ])
                    for lab, fmt in [
                        ('State Loop', 'movie.add_state_loop(%d, %d'),
                        ('State Sweep', 'movie.add_state_sweep(%d, %d'),
                    ]
                ]),
                ('command', 'Update Last Program', self.mvprg),
                ('command', 'Remove Last Program', self.mvprg_remove_last),
                ('separator',),
                ('command', 'Reset', 'mset;rewind'),
                ('separator',),
                ('menu', 'Frame Rate', [
                    ('radio', '30 FPS',  'movie_fps', 30.),
                    ('radio', '15 FPS',  'movie_fps', 15.),
                    ('radio', '5 FPS',   'movie_fps',  5.),
                    ('radio', '1 FPS',   'movie_fps',  1.),
                    ('radio', '0.3 FPS', 'movie_fps',  .3),
                    ('separator',),
                    ('check', 'Show Frame Rate', 'show_frame_rate'),
                    ('command', 'Reset Meter', cmd.meter_reset),
                ]),
                ('separator',),
                ('check', 'Auto Interpolate', 'movie_auto_interpolate'),
                ('check', 'Show Panel', 'movie_panel'),
                ('check', 'Loop Frames', 'movie_loop'),
                ('check', 'Draw Frames', 'draw_frames'),
                ('check', 'Ray Trace Frames', 'ray_trace_frames'),
                ('check', 'Cache Frame Images', 'cache_frames'),
                ('command', 'Clear Image Cache', cmd.mclear),
                ('separator',),
                ('check', 'Static Singletons', 'static_singletons', 1),
                ('check', 'Show All States', 'all_states', 1),
            ]),
            ('menu', 'Display', [
                ('check', 'Sequence', 'seq_view', 1),
                ('menu', 'Sequence Mode', [
                    ('radio', lab, 'seq_view_format', val)
                    for lab, val in [
                        ('Residue Codes', 0),
                        ('Residue Names', 1),
                        ('Chain Identifiers', 3),
                        ('Atom Names', 2),
                        ('States', 4),
                    ]
                ] + [
                    ('separator',),
                ] + [
                    ('radio', lab, 'seq_view_label_mode', val)
                    for lab, val in [
                        ('All Residue Numbers', 2),
                        ('Top Sequence Only', 1),
                        ('Object Names Only', 0),
                        ('No Labels', 3),
                    ]
                ] + [
                    ('separator',),
                ] + [
                    ('radio', lab, 'seq_view_gap_mode', val)
                    for lab, val in [
                        ('No Gaps', 0),
                        ('All Gaps', 1),
                        ('Single Gap', 2),
                    ]
                ]),
                ('separator',),
                ('check', 'Internal GUI', 'internal_gui', 1),
                ('check', 'Internal Prompt', 'internal_prompt', 1),
                ('menu', 'Internal Feedback', [
                    ('radio', str(val), 'internal_feedback', val)
                    for val in [0, 1, 3, 5]
                ]),
                ('menu', 'Overlay', [
                    ('radio', str(val), 'overlay', val)
                    for val in [0, 1, 3, 5]
                ]),
                ('separator',),
                ('check', 'Stereo', 'stereo', 1),
                ('menu', 'Stereo Mode', [
                    ('command', 'Anaglyph Stereo', 'stereo anaglyph'),
                    ('command', 'Cross-Eye Stereo', 'stereo crosseye'),
                    ('command', 'Wall-Eye Stereo', 'stereo walleye'),
                    ('command', 'Quad-Buffered Stereo', 'stereo quadbuffer'),
                    ('command', 'Zalman Stereo', 'stereo byrow'),
                    ('command', 'OpenVR', 'stereo openvr'),
                    ('separator',),
                    ('command', 'Swap Sides', 'stereo swap'),
                    ('separator',),
                    ('command', 'Chromadepth', 'stereo chromadepth'),
                    ('command', 'off', 'stereo off'),
                ]),
                ('separator',),
                ('menu', 'Zoom', [
                    ('command', str(i) + ' Angstrom Sphere', lambda i=i: cmd.zoom('center', i, animate=-1))
                    for i in [4, 6, 8, 12, 20]
                ] + [
                    ('command', 'All', 'zoom animate=-1'),
                    ('command', 'Complete', 'zoom animate=-1, complete=1'),
                ]),
                ('menu', 'Clip', [
                    ('command', 'Nothing', 'clip atoms, 5, all'),
                ] + [
                    ('command', str(i) + ' Angstrom Slab', lambda i=i: cmd.clip('slab', i))
                    for i in [8, 12, 16, 20, 30]
                ]),
                ('separator',),
                ('menu', 'Background', [
                    ('check', 'Opaque', 'opaque_background', 1),
                    ('check', 'Alpha Checker', 'show_alpha_checker', 1),
                    ('separator',),
                ] + [
                    ('radio', lab, 'bg_rgb', val)
                    for lab, val in [
                        ('White', 0),           # white
                        ('Light Grey', 134),    # grey80
                        ('Grey', 104),          # grey50
                        ('Black', 1),           # black
                    ]
                ]),
                ('menu', 'Color Space', [
                    ('command', 'CMYK (for publications)',  'space cmyk'),
                    ('command', 'PyMOL (for video + web)',  'space pymol'),
                    ('command', 'RGB (default)',            'space rgb'),
                ]),
                ('menu', 'Quality', [
                    ('command', 'Maximum Performance',      'util.performance(100)'),
                    ('command', 'Reasonable Performance',   'util.performance(66)'),
                    ('command', 'Reasonable Quality',       'util.performance(33)'),
                    ('command', 'Maximum Quality',          'util.performance(0)'),
                ]),
                ('menu', 'Grid', [
                    ('radio', lab, 'grid_mode', val)
                    for lab, val in [
                        ('By Object', 1),
                        ('By State', 2),
                        ('By Object-State', 3),
                        ('Disable', 0),
                    ]
                ]),
                ('separator',),
                ('check', 'Orthoscopic View', 'orthoscopic', 1),
                ('check', 'Show Valences', 'valence', 1),
                ('check', 'Smooth Lines', 'line_smooth', 1),
                ('check', 'Depth Cue (Fogging)', 'depth_cue', 1),
                ('check', 'Two Sided Lighting', 'two_sided_lighting', 1),
                ('check', 'Specular Reflections', 'specular', 1.0), # offvalue=0.0,
                ('check', 'Animation', 'animation', 1),
                ('check', 'Roving Detail', 'roving_detail', 1),
            ]),
            ('menu', 'Setting', [
                ('command', 'Edit All...', self.settings_edit_all_dialog),
                ('command', 'Colors...', self.edit_colors_dialog),
                ('separator',),
                ('menu', 'Label', [
                    ('menu', 'Size', [
                        ('radio', str(val) + ' Point', 'label_size', val)
                        for val in [10, 14, 18, 24, 36, 48, 72]
                    ] + [
                        ('separator',),
                    ] + [
                        ('radio', str(val) + ' Angstrom', 'label_size', -val)
                        for val in [0.3, 0.5, 1, 2, 4]
                    ]),
                    ('menu', 'Font', [
                        ('radio', label, 'label_font_id', val)
                        for label, val in [
                            ('Sans', 5),
                            ('Sans Oblique', 6),
                            ('Sans Bold', 7),
                            ('Sans Bold Oblique', 8),
                            ('Serif', 9),
                            ('Serif Oblique',17),
                            ('Serif Bold', 10),
                            ('Serif Bold Oblique', 18),
                            ('Mono', 11),
                            ('Mono Oblique', 12),
                            ('Mono Bold', 13),
                            ('Mono Bold Oblique', 14),
                            ('Gentium Roman', 15),
                            ('Gentium Italic', 16),
                        ]
                    ]),
                    ('menu', 'Color', [
                        ('radio', lab, 'label_color', val)
                        for lab, val in [
                            ('Front', -6),
                            ('Back', -7),
                        ]
                    ]),
                    ('check', 'Show Connectors', 'label_connector'),
                    ('menu', 'Background Color', [
                        ('radio', lab, 'label_bg_color', val)
                        for lab, val in [
                            ('None', -1),
                            ('Back', -7),
                            ('Front', -6),
                        ]
                    ]),
                ]),
                ('menu', 'Lines & Sticks', [
                    ('check', 'Ball and Stick', 'stick_ball', 1),
                    ('menu', 'Ball and Stick Ratio', [
                        ('radio', lab, 'stick_ball_ratio', val)
                        for lab, val in [
                            ('1.0', 1.),
                            ('1.5', 1.5),
                            ('VDW', -1.),
                        ]
                    ]),
                    ('separator',),
                    ('menu', 'Zero Order Bonds', [
                        ('radio', lab, 'valence_zero_mode', val)
                        for lab, val in [
                            ('Hide', 0),
                            ('Dashed', 1),
                            ('Solid', 2),
                        ]
                    ]),
                    ('menu', 'Zero Order Stick Scale', [
                        ('radio', str(val), 'valence_zero_scale', val)
                        for val in [0.1, 0.2, 0.3, 1.0]
                    ]),
                    ('separator',),
                    ('menu', 'Stick Radius', [
                        ('radio', str(val), 'stick_radius', val)
                        for val in [.1, .2, .25]
                    ]),
                    ('menu', 'Stick Hydrogen Scale', [
                        ('radio', str(val), 'stick_h_scale', val)
                        for val in [.4, 1.]
                    ]),
                    ('separator',),
                    ('menu', 'Line Width', [
                        ('radio', str(val), 'line_width', val)
                        for val in [1.0, 1.49, 3.0]
                    ]),
                    ('check', 'Lines As Cylinders', 'line_as_cylinders', 1),
                ]),
                ('menu', 'Cartoon', [
                    ('menu', 'Rings and Bases', [
                        ('radio', 'Filled Rings (Round Edges)', 'cartoon_ring_mode', 1),
                        ('radio', 'Filled Rings (Flat Edges)', 'cartoon_ring_mode', 2),
                        ('radio', 'Filled Rings (with Border)', 'cartoon_ring_mode', 3),
                        ('radio', 'Spheres', 'cartoon_ring_mode', 4),
                        ('radio', 'Base Ladders', 'cartoon_ring_mode', 0),
                        ('separator',),
                        ('radio', 'Bases and Sugars', 'cartoon_ring_finder', 1),
                        ('radio', 'Bases Only', 'cartoon_ring_finder', 2),
                        ('radio', 'Non-protein Rings', 'cartoon_ring_finder', 3),
                        ('radio', 'All Rings', 'cartoon_ring_finder', 4),
                        ('separator',),
                        ('radio', 'Transparent Rings', 'cartoon_ring_transparency', 0.5),
                        ('radio', 'Default', 'cartoon_ring_transparency', -1),
                    ]),
                    ('check', 'Side Chain Helper', 'cartoon_side_chain_helper', 1),
                    ('check', 'Round Helices', 'cartoon_round_helices', 1),
                    ('check', 'Fancy Helices', 'cartoon_fancy_helices', 1),
                    ('check', 'Cylindrical Helices', 'cartoon_cylindrical_helices', 1),
                    ('check', 'Flat Sheets', 'cartoon_flat_sheets', 1),
                    ('check', 'Fancy Sheets', 'cartoon_fancy_sheets', 1),
                    ('check', 'Smooth Loops', 'cartoon_smooth_loops', 1),
                    ('check', 'Discrete Colors', 'cartoon_discrete_colors', 1),
                    ('check', 'Highlight Color', 'cartoon_highlight_color', 104, -1),
                    ('menu', 'Sampling', [
                        ('radio', 'Atom count dependent', 'cartoon_sampling', -1),
                    ] + [
                        ('radio', str(val), 'cartoon_sampling', val)
                        for val in [2, 7, 14]
                    ]),
                    ('menu', 'Gap Cutoff', [
                        ('radio', str(val), 'cartoon_gap_cutoff', val)
                        for val in [0, 5, 10, 20]
                    ]),
                ]),
                ('menu', 'Ribbon', [
                    ('check', 'Side Chain Helper', 'ribbon_side_chain_helper', 1),
                    ('check', 'Trace Atoms', 'ribbon_trace_atoms', 1),
                    ('separator',),
                    ('radio', 'As Lines', 'ribbon_as_cylinders', 0),
                    ('radio', 'As Cylinders', 'ribbon_as_cylinders', 1),
                    ('menu', 'Cylinder Radius', [
                        ('radio', 'Match Line Width', 'ribbon_radius', 0.),
                    ] + [
                        ('radio', '%.1f Angstrom' % val, 'ribbon_radius', val)
                        for val in [.2, .5, 1.]
                    ]),
                ]),
                ('menu', 'Surface', [
                    ('menu', 'Color', [
                        ('radio', 'White', 'surface_color', 0), # white
                        ('radio', 'Light Gray', 'surface_color', 4236), # gray80
                        ('radio', 'Gray', 'surface_color', 25), # gray
                        ('radio', 'Default (Atomic)', 'surface_color', -1),
                    ]),
                    ('radio', 'Dot', 'surface_type', 1),
                    ('radio', 'Wireframe', 'surface_type', 2),
                    ('radio', 'Solid', 'surface_type', 0),
                    ('separator',),
                    ('radio', 'Cavities and Pockets Only', 'surface_cavity_mode', 1),
                    ('radio', 'Cavities and Pockets (Culled)', 'surface_cavity_mode', 2),
                    ('menu', 'Cavity Detection Radius', [
                        ('radio', str(val) + ' Angstrom', 'surface_cavity_radius', val)
                        for val in [7]
                    ] + [
                        ('radio', str(val) + ' Solvent Radii', 'surface_cavity_radius', -val)
                        for val in [3, 4, 5, 6, 8, 10, 20]
                    ]),
                    ('menu', 'Cavity Detection Cutoff', [
                        ('radio', str(val) + ' Solvent Radii', 'surface_cavity_cutoff', -val)
                        for val in [1, 2, 3, 4, 5]
                    ]),
                    ('radio', 'Exterior (Normal)', 'surface_cavity_mode', 0),
                    ('separator',),
                    ('check', 'Solvent Accessible', 'surface_solvent', 1),
                    ('separator',),
                    ('check', 'Smooth Edges', 'surface_smooth_edges', 1),
                    ('check', 'Edge Proximity', 'surface_proximity', 1),
                    ('separator',),
                    ('radio', 'Ignore None', 'surface_mode', 1),
                    ('radio', 'Ignore HETATMs', 'surface_mode', 0),
                    ('radio', 'Ignore Hydrogens', 'surface_mode', 2),
                    ('radio', 'Ignore Unsurfaced', 'surface_mode', 3),
                ]),
                ('menu', 'Volume', [
                    ('check', 'Pre-integrated Rendering', 'volume_mode'),
                    ('menu', 'Number of Layers', [
                        ('radio', '%.0f' % val, 'volume_layers', val)
                        for val in (100., 256., 500., 1000.)
                    ]),
                ]),
                ('menu', 'Transparency', [
                    ('menu', 'Surface', transparency_menu('transparency')),
                    ('menu', 'Sphere',  transparency_menu('sphere_transparency')),
                    ('menu', 'Cartoon', transparency_menu('cartoon_transparency')),
                    ('menu', 'Stick',   transparency_menu('stick_transparency')),
                    ('separator',),
                ] + [
                    ('command', lab, lambda v=val: (
                        cmd.set('transparency_mode', v[0], quiet=0),
                        cmd.set('backface_cull', v[1], quiet=0),
                        cmd.set('two_sided_lighting', v[2], quiet=0)))
                    for lab, val in [
                        ('Uni-Layer',     (2, 1, 0)),
                        ('Multi-Layer',   (1, 0, 1)),
                        ('Multi-Layer (Real-time OIT)', (3, 0, -1)),
                        ('Fast and Ugly', (0, 1, 0)),
                    ]
                ] + [
                    ('separator',),
                    ('check', 'Angle-dependent', 'ray_transparency_oblique', 1.0), # offvalue=0.0,
                ]),
                ('menu', 'Rendering', [
                    ('check', 'OpenGL 2.0 Shaders', 'use_shaders', 1),
                    ('separator',),
                    ('check', 'Antialias (Ray Tracing)', 'antialias', 1),
                    ('menu', 'Antialias (Real Time)', [
                        ('radio', lab, 'antialias_shader', val)
                        for lab, val in [('off', 0), ('FXAA', 1), ('SMAA', 2)]
                    ]),
                    ('separator',),
                    ('command', 'Modernize', lambda: cmd.util.modernize_rendering(1, cmd)),
                    ('separator',),
                    ('menu', 'Shadows', [
                        ('command', val.title(), lambda v=val: cmd.util.ray_shadows(v))
                        for val in ['none', 'light', 'medium', 'heavy', 'black']
                    ] + [
                        ('separator',),
                    ] + [
                        ('command', val.title(), lambda v=val: cmd.util.ray_shadows(v))
                        for val in ['matte', 'soft', 'occlusion', 'occlusion2']
                    ]),
                    ('menu', 'Texture', [
                        ('radio', 'None', 'ray_texture', 0),
                        ('radio', 'Matte 1', 'ray_texture', 1),
                        ('radio', 'Matte 2', 'ray_texture', 4),
                        ('radio', 'Swirl 1', 'ray_texture', 2),
                        ('radio', 'Swirl 2', 'ray_texture', 3),
                        ('radio', 'Fiber', 'ray_texture', 5),
                    ]),
                    ('menu', 'Interior Texture', [
                        ('radio', 'None', 'ray_interior_texture', 0),
                        ('radio', 'Matte 1', 'ray_interior_texture', 1),
                        ('radio', 'Matte 2', 'ray_interior_texture', 4),
                        ('radio', 'Swirl 1', 'ray_interior_texture', 2),
                        ('radio', 'Swirl 2', 'ray_interior_texture', 3),
                        ('radio', 'Fiber', 'ray_interior_texture', 5),
                    ]),
                    ('menu', 'Memory', [
                        ('radio', lab, 'hash_max', val)
                        for lab, val in [
                            ('Use Less (slower)', 70),
                            ('Use Standard Amount', 100),
                            ('Use More (faster)', 170),
                            ('Use Even More', 230),
                            ('Use Most', 300),
                        ]
                    ]),
                    ('separator',),
                    ('check', 'Cull Backfaces', 'backface_cull', 1),
                    ('check', 'Opaque Interiors', 'ray_interior_color', 74, -1),
                ]),
                ('separator',),
                ('menu', 'PDB File Loading', [
                    ('check', 'Ignore PDB Segment Identifier', 'ignore_pdb_segi', 1),
                ]),
                ('menu', 'mmCIF File Loading', [
                    ('check', 'Use "auth" Identifiers', 'cif_use_auth', 1),
                    ('check', 'Load Assembly (Biological Unit)', 'assembly', "1", ""),
                    ('check', 'Bonding by "Chemical Component Dictionary"', 'connect_mode', 4, 0),
                ]),
                ('menu', 'Map File Loading', [
                    ('check', 'Normalize CCP4 Maps', 'normalize_ccp4_maps', 1),
                    ('check', 'Normalize O Maps', 'normalize_o_maps', 1),
                ]),
                ('separator',),
                ('menu', 'Auto-Show ...', [
                    ('check', 'Cartoon/Sticks/Spheres by Classification', 'auto_show_classified', -1, 0),
                    ('separator',),
                    ('check', 'Auto-Show Lines', 'auto_show_lines', 1),
                    ('check', 'Auto-Show Spheres', 'auto_show_spheres', 1),
                    ('check', 'Auto-Show Nonbonded', 'auto_show_nonbonded', 1),
                    ('separator',),
                    ('check', 'Auto-Show New Selections', 'auto_show_selections', 1),
                    ('check', 'Auto-Hide Selections', 'auto_hide_selections', 1),
                ]),
                ('check', 'Auto-Zoom New Objects', 'auto_zoom', 1),
                ('check', 'Auto-Remove Hydrogens', 'auto_remove_hydrogens', 1),
                ('separator',),
                ('check', 'Show Text (Esc)', 'text'),
                ('check', 'Overlay Text', 'overlay'),
            ]),
            ('menu', 'Scene', [
                ('command', 'Next [PgDn]', lambda: cmd.scene('', 'next')),
                ('command', 'Previous [PgUp]', lambda: cmd.scene('', 'previous')),
                ('separator',),
                ('command', 'Append', 'scene new, store'),
                ('menu', 'Append...', [
                    ('command', 'Camera', 'scene new, store, color=0, rep=0'),
                    ('command', 'Color', 'scene new, store, view=0, rep=0'),
                    ('command', 'Reps', 'scene new, store, view=0, color=0'),
                    ('command', 'Reps + Color', 'scene new, store, view=0'),
                ]),
                ('command', 'Insert Before', lambda: cmd.scene('', 'insert_before')),
                ('command', 'Insert After', lambda: cmd.scene('','insert_after')),
                ('command', 'Update', lambda: cmd.scene('auto','update')),
                ('separator',),
                ('command', 'Delete', lambda: cmd.scene('auto','clear')),
                ('separator',),
                ('menu', 'Recall', F_scene_menu('recall')),
                ('menu', 'Store', F_scene_menu('store')),
                ('menu', 'Clear', F_scene_menu('clear')),
                ('separator',),
                ('check', 'Buttons', 'scene_buttons', 1),
                ('menu', 'Cache', [
                    ('command', 'Enable', lambda: cmd.cache("enable")),
                    ('command', 'Optimize', lambda: cmd.cache("optimize")),
                    ('command', 'Read Only', lambda: cmd.cache("read_only")),
                    ('command', 'Disable', lambda: cmd.cache("disable")),
                ]),
            ]),
            ('menu', 'Mouse', [
                ('menu', 'Selection Mode', [
                    ('radio', lab, 'mouse_selection_mode', val)
                    for lab, val in [
                        ('Atoms', 0),
                        ('Residues', 1),
                        ('Chains', 2),
                        ('Segments', 3),
                        ('Objects', 4),
                        ('Molecules', 5),
                        ('C-alphas', 6),
                    ]
                ]),
                ('separator',),
                ('command', '3 Button Motions', lambda: cmd.config_mouse('three_button_motions')),
                ('command', '3 Button Editing', lambda: cmd.config_mouse('three_button_editing')),
                ('command', '3 Button Viewing', lambda: cmd.mouse('three_button_viewing')),
                ('command', '3 Button Lights', lambda: cmd.mouse('three_button_lights')),
                ('command', '3 Button All Modes', lambda: cmd.config_mouse('three_button_all_modes')),
                ('command', '2 Button Editing', lambda: cmd.config_mouse('two_button_editing')),
                ('command', '2 Button Viewing', lambda: cmd.config_mouse('two_button')),
                ('command', '1 Button Viewing Mode', lambda: cmd.mouse('one_button_viewing')),
                ('command', 'Emulate Maestro', lambda: cmd.mouse('three_button_maestro')),
                ('separator',),
                ('check', 'Virtual Trackball', 'virtual_trackball'),
                ('check', 'Show Mouse Grid', 'mouse_grid'),
                ('check', 'Roving Origin', 'roving_origin'),
            ]),
            ('menu', 'Wizard', [
                ('command', 'Appearance', 'wizard appearance'),
                ('command', 'Measurement', 'wizard measurement'),
                ('menu', 'Mutagenesis', [
                    ('command', 'Protein', 'wizard mutagenesis'),
                    ('command', 'Nucleic Acids', 'wizard nucmutagenesis'),
                ]),
                ('command', 'Pair Fitting', 'wizard pair_fit'),
                ('separator',),
                ('command', 'Density', 'wizard density'),
                ('command', 'Filter', 'wizard filter'),
                ('command', 'Sculpting', 'wizard sculpting'),
                ('separator',),
                ('command', 'Label', 'wizard label'),
                ('command', 'Charge', 'wizard charge'),
                ('separator',),
                ('menu', 'Demo', [
                    ('command', 'Representations', lambda: cmd.wizard('demo', 'reps')),
                    ('command', 'Cartoon Ribbons', lambda: cmd.wizard('demo', 'cartoon')),
                    ('command', 'Roving Detail', lambda: cmd.wizard('demo', 'roving')),
                    ('command', 'Roving Density', lambda: cmd.wizard('demo', 'roving_density')),
                    ('command', 'Transparency', lambda: cmd.wizard('demo', 'trans')),
                    ('command', 'Ray Tracing', lambda: cmd.wizard('demo', 'ray')),
                    ('command', 'Sculpting', lambda: cmd.wizard('demo', 'sculpt')),
                    ('command', 'Scripted Animation', lambda: cmd.wizard('demo', 'anime')),
                    ('command', 'Electrostatics', lambda: cmd.wizard('demo', 'elec')),
                    ('command', 'Compiled Graphics Objects', lambda: cmd.wizard('demo', 'cgo')),
                    ('command', 'Molscript/Raster3D Input', lambda: cmd.wizard('demo', 'raster3d')),
                    ('separator',),
                    ('command', 'End Demonstration', lambda: cmd.replace_wizard('demo', 'finish')),
                ]),
            ]),
            ('menu', 'Plugin', []),
            ('menu', 'Help', [
                ('command', 'PyMOL Command Reference', lambda: webbrowser.open('http://pymol.org/pymol-command-ref.html')),
                ('separator',),
                ('command', 'Online Documentation', lambda: webbrowser.open("http://pymol.org/d/")),
                ('menu', 'Topics', [
                    ('command', 'Introductory Screencasts', lambda: webbrowser.open("http://pymol.org/d/media:intro")),
                    ('command', 'Core Commands', lambda: webbrowser.open("http://pymol.org/d/command:core_set")),
                    ('command', 'Settings', lambda: webbrowser.open("http://pymol.org/d/setting")),
                    ('command', 'Atom Selections', lambda: webbrowser.open("http://pymol.org/d/selection")),
                    ('command', 'Commands', lambda: webbrowser.open("http://pymol.org/d/command")),
                    ('command', 'Launching', lambda: webbrowser.open("http://pymol.org/d/launch")),
                    ('command', 'Concepts', lambda: webbrowser.open("http://pymol.org/d/concept")),
                    ('command', 'A.P.I. Methods', lambda: webbrowser.open("http://pymol.org/d/api")),
                ]),
                ('separator',),
                ('command', 'PyMOL Community Wiki', lambda: webbrowser.open("http://www.pymolwiki.org")),
                ('command', 'PyMOL Mailing List', lambda: webbrowser.open("https://lists.sourceforge.net/lists/listinfo/pymol-users")),
                ('command', 'PyMOL Home Page', lambda: webbrowser.open("http://www.pymol.org")),
                ('separator',),
                ('command', 'About PyMOL', self.show_about),
                ('command', 'Sponsorship Information', lambda: webbrowser.open("http://pymol.org/funding.html")),
                ('command', 'How to Cite PyMOL', lambda: webbrowser.open("http://pymol.org/citing")),
            ]),
        ]

    #################
    # command history
    #################

    def _setup_history(self):
        self.history = ['']
        self.history_cur = 0

    def complete(self, *args):
        st = self.cmd._parser.complete(self.command_get())
        if st:
            self.command_set(st)
            self.command_set_cursor(len(st))
        return 'break'

    def doTypedCommand(self, cmmd):
        '''Run command from command line'''
        if len(self.history) < 2 or self.history[1] != cmmd:
            self.history[0] = cmmd
            self.history.insert(0, '') # always leave blank at 0
            if len(self.history) > 255:
                self.history.pop()
        self.history_cur = 0
        self.cmd.do(cmmd)

    def back_search(self, set0=False):
        '''Command line history back search'''
        if not self.history_cur or set0:
            self.history[0] = self.command_get()
        for i in range(self.history_cur + 1, len(self.history)):
            if self.history[i].startswith(self.history[0]):
                self._jump_history(i)
                break

    def back(self):
        '''Command line history step back'''
        if not self.history_cur:
            self.history[0] = self.command_get()
        self._jump_history(self.history_cur + 1)

    def forward(self):
        '''Command line history step forward'''
        if not self.history_cur:
            return
        self._jump_history(self.history_cur - 1)

    def _jump_history(self, i):
        self.history_cur = min(i, len(self.history) - 1)
        self.command_set(self.history[self.history_cur])
        l = len(self.history[self.history_cur])
        self.command_set_cursor(l)

    #######################
    # Movie Program helpers
    #######################

    movie_start = 0
    movie_command = None

    def mvprg_remove_last(self):
        '''
        Remove frames [N..end] where N is the beginning of the last added
        program.
        '''
        if self.movie_start > 0:
            self.cmd.mdelete(-1, self.movie_start)

    def mvprg(self, command=None):
        '''
        Movie program helper function with "Update" and "Remove" support
        for the last added program
        '''
        if command is not None:
            self.movie_start = self.cmd.get_movie_length() + 1
            self.movie_command = command % self.movie_start
        elif self.movie_command is None:
            return

        self.cmd.do(self.movie_command)

    #######################
    # Recent files history
    #######################

    _recent_filenames_db = None

    def _recent_filenames_lazy_init(self):
        '''Lazy initialization of the recent files database. Returns True
        on success or if the database has already been initialized.
        '''
        if self._recent_filenames_db is not None:
            if self._recent_filenames_db is False:
                return False
            return True

        d = os.path.expanduser('~/.pymol')
        f = os.path.join(d, 'recent.db')

        try:
            os.makedirs(d)
        except OSError:
            pass

        try:
            import sqlite3
            db = sqlite3.connect(f)
            db.cursor().execute('CREATE TABLE IF NOT EXISTS '
                    'recent (filename text unique, timestamp integer)')
            self._recent_filenames_db = db
        except BaseException as e:
            print(' Warning: failed to connect to recent DB:', e)
            self._recent_filenames_db = False
            return False

        return True

    @property
    def recent_filenames(self):
        '''List of recently opened files'''
        if not self._recent_filenames_lazy_init():
            return []
        c = self._recent_filenames_db.cursor()
        return [row[0] for row in c.execute(
            'SELECT filename FROM recent ORDER BY timestamp DESC')]

    def recent_filenames_add(self, filename):
        '''Register a file name as a "recently opened file"'''
        if not self._recent_filenames_lazy_init():
            return

        c = self._recent_filenames_db.cursor()
        c.execute("REPLACE INTO recent VALUES (?, datetime('now'))", (filename,))

        # keep number reasonable (this could be a one-liner if sqlite is
        # compiled with SQLITE_ENABLE_UPDATE_DELETE_LIMIT)
        if c.execute('SELECT COUNT(*) FROM recent').fetchone()[0] > 20:
            c.execute('DELETE FROM recent WHERE timestamp < ?',
                    c.execute('SELECT timestamp from recent '
                        'ORDER BY timestamp DESC '
                        'LIMIT 1 OFFSET 15').fetchone())

        self._recent_filenames_db.commit()
