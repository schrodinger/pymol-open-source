files_names = [
      'champ.c',
      'champ_module.c',
      'chiral.c',
      'err2.c',
      'feedback2.c',
      'list.c',
      'os_memory.c',
      'sort.c',
      'strblock.c',
      'vla.c',
]
champ_src_files = []

foreach file_name : files_names
  champ_src_files += files(file_name)
endforeach

py.extension_module(
  '_champ',
  champ_src_files,
  include_directories : include_directories('.'),
  install: true,
  subdir: 'chempy/champ'
)

