'''
Update PyMOL's copy of the molfile plugins
'''

import os

molfile_src_path = "/tmp/plugins/molfile_plugin/src"

if not os.path.isdir(molfile_src_path):
    raise SystemExit('noch such dir: ' + molfile_src_path)

# remove existing generated files
os.system("/bin/rm src/*")

src_list = []

blacklist = [
    'babelplugin',      # requires openbabel
    'cpmdlogplugin',    #
    'cpmdplugin',       #
    'gaussianplugin',   #
    'hoomdplugin',      # requires expat
    'lammpsplugin',     # requires gz
    'webpdbplugin',     # tcl dependent
    'tngplugin',        # requires Gromacs TNG library
    'dmsplugin',        # requires sqlite3
]

for in_file in sorted(os.listdir(molfile_src_path)):
    pref, ext = os.path.splitext(in_file)

    if not pref.endswith('plugin'):
        continue

    if pref in blacklist:
        continue

    if ext in ['.C', '.cxx']:
        ext = '.cpp'
    elif ext not in ['.c']:
        continue

    print("processing: " + pref)

    in_file = os.path.join(molfile_src_path, in_file)
    out_file = os.path.join('src', pref + ext)

    input = open(in_file, 'rb').readlines()

    src_list.append(pref)

    with open(out_file,'wb') as g:
        g.write(b"/* MACHINE GENERATED FILE, DO NOT EDIT! */\n\n")
        g.write(b"#define VMDPLUGIN molfile_%s\n" % pref.encode())
        g.write(b"#define STATIC_PLUGIN 1\n\n")
        for i, line in enumerate(input, 1):
            # no including of hash.c, inthash.c like header files
            if line.startswith(b'#define VMDPLUGIN_STATIC'):
                continue

            # included, don't compile separatly
            line = line.replace(b'"ply.c"', b'"ply_c.h"')

            g.write(line)

with open("src/PlugIOManagerInit.c", 'w') as g:
    g.write("/* MACHINE GENERATED FILE, DO NOT EDIT! */\n\n")
    g.write('#include "vmdplugin.h"\n\n')
    g.write('struct PyMOLGlobals;\n');
    g.write('/* prototypes */\n')
    for pref in src_list:
        g.write("int molfile_%s_init(void);\n"%pref)
        g.write("int molfile_%s_register(void *,vmdplugin_register_cb);\n"%pref)
        g.write("int molfile_%s_fini(void);\n"%pref)
    g.write('''

    int PlugIOManagerRegister(struct PyMOLGlobals *G, vmdplugin_t *);

    int PlugIOManagerInitAll(struct PyMOLGlobals *G);

    int PlugIOManagerInitAll(struct PyMOLGlobals *G)
    {
       int ok=1;
''')
    for pref in src_list:
        g.write("if(ok) ok = ok && (molfile_%s_init() == VMDPLUGIN_SUCCESS);\n"%pref)
    g.write('''
       if(ok) {
''')
    for pref in src_list:
        g.write("if(ok) ok = ok && (molfile_%s_register(G,(vmdplugin_register_cb)PlugIOManagerRegister) == VMDPLUGIN_SUCCESS);\n"%pref)
    g.write('''
       }
       return ok;
    }
    ''')
    g.write('''

    int PlugIOManagerFreeAll(void);
    int PlugIOManagerFreeAll(void)
    {
       int ok=1;
''')
    for pref in src_list:
        g.write("if(ok) ok = ok && (molfile_%s_fini() == VMDPLUGIN_SUCCESS);\n"%pref)

    g.write('''
       return ok;
    }\n''')


if True:
    os.system("/bin/cp %s/*.h* src/"%molfile_src_path)
    os.system("/bin/cp %s/hash.c src/"%molfile_src_path)
    os.system("/bin/cp %s/inthash.c src/" % molfile_src_path)
    os.system("/bin/cp %s/ply.c src/ply_c.h" % molfile_src_path) # included
    os.system("/bin/cp %s/../../include/*.h ../include/" % molfile_src_path)
    os.system("/bin/cp %s/../LICENSE ./" % molfile_src_path)

    os.system("patch -p5 -i post.patch")

    os.system("/bin/chmod -x src/*")
