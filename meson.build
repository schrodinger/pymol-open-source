project(
  'pymol',
  'c', 'cpp',
  version: run_command(
    ['modules/version.py'],
    check: true
  ).stdout().strip(),
  meson_version: '>= 1.4.0',
  default_options: [
    'cpp_std=c++17',
  ],
)

py = import('python').find_installation(pure: false)

dependencies = [
  dependency('libpng', fallback: 'libpng'),
  dependency('freetype2', fallback: 'freetype'),
  dependency('gl', 'opengl'),
  dependency('glew', fallback: 'glew'),
  dependency('glm', fallback: 'glm'),
]

run_command(py, files('create_shadertext.py'))
src_files = files('build/ShaderText.cpp')
include_dirs = [
  include_directories('include'),
  include_directories('build'),
]
incdir_numpy = run_command(py,
  [
    '-c',
    'import os, numpy; print(os.path.relpath(numpy.get_include()))'
  ],
  check: true
).stdout().strip()
include_dirs += include_directories(incdir_numpy)

macros = [
  '-D_PYMOL_LIBPNG',
  '-D_PYMOL_FREETYPE',
  '-D_PYMOL_NUMPY',
]
flags = host_machine.system() == 'windows' ? ['-WMP'] : [
  '-Werror=return-type',
  '-Wunused-variable',
  '-Wno-switch',
  '-Wno-narrowing',
  '-Wno-char-subscripts',
]


#if get_option('use-openmp')
#  dependencies += dependency('openmp')
#  # shouldn't it be _PYMOL_OPENMP?
#  macros += '-DPYMOL_OPENMP'
#endif

if not get_option('no-vmd-plugins')
  include_dirs += include_directories('contrib/uiuc/plugins/include')
  src_dirs += 'contrib/uiuc/plugins/molfile_plugin/src'
  dependencies += dependency('netcdf')
  macros += '-D_PYMOL_VMD_PLUGINS'
endif

if not get_option('no-libxml')
  dependencies += dependency('libxml-2.0')
  macros += '-D_HAVE_LIBXML'
endif

msg_dep = get_option('use-msgpackc')
if msg_dep != 'guess'
  dependencies += dependency('msgpackc')
else
  macros += '-D_PYMOL_NO_MSGPACKC'
endif

if not get_option('no-glut')
  macros += '-D_PYMOL_NO_MAIN'
endif

if get_option('testing')
  subdir('layerCTest')
endif

if get_option('openvr')
  dependencies += dependency('openvr')
  macros += '-D_PYMOL_OPENVR'
endif

# macos
if host_machine.system() == 'darwin'
  dependencies += [
    dependency('glut'),
  ]
  # shouldn't it be _PYMOL_CURVE_VALIDATE?
  macros += '-DPYMOL_CURVE_VALIDATE'

  if get_option('no-osx-frameworks')
      macros += '-D_PYMOL_OSX'
  endif
endif

if host_machine.system() == 'linux'
  dependencies += [
    dependency('glut'),
  ]
endif


subdir('ov/src')
subdir('layer0')
subdir('layer1')
subdir('layer2')
subdir('layer3')
subdir('layer4')
subdir('layer5')

py.extension_module(
  '_cmd',
  src_files,
  include_directories : include_dirs,
  cpp_args: macros + flags,
  dependencies: dependencies,
  install: true,
  subdir: 'pymol'
)
subdir('contrib/champ')

subdir('modules/chempy')
subdir('modules/pmg_qt')
subdir('modules/pmg_tk')
subdir('modules/pymol')
subdir('modules/pymol2')

