if 1:
    print "DANGER DO NOT RUN UNTIL changes to Gromacs.h and gromacsplugin.ccp are handled"
else :

    import re
    import string
    import os
    from glob import glob

    molfile_src_path = "/home/warren/06duo/software/vmd/plugins/molfile_plugin/src"

    src_list=[
        'avsplugin',
    #    'babelplugin',
        'bgfplugin',
        'binposplugin',
        'biomoccaplugin',
        'brixplugin',
        'carplugin',
        'ccp4plugin',
        'corplugin',
        'cpmdplugin',
        'crdplugin',
        'cubeplugin',
        'dcdplugin',
        'dlpolyplugin',
        'dsn6plugin',
        'dxplugin',
        'edmplugin',
        'fs4plugin',
        'gamessplugin',
        'graspplugin',
        'grdplugin',
        'gridplugin',
        'gromacsplugin',
    ##    'lammpsplugin',
        'mapplugin',
        'mdfplugin',
        'mol2plugin',
        'moldenplugin',
        'msmsplugin',
        'namdbinplugin',
        'parm7plugin',
        'parmplugin',
        'pdbplugin',
        'phiplugin',
        'pltplugin',
        'pqrplugin',
        'psfplugin',
        'raster3dplugin',
        'rst7plugin',
        'situsplugin',
        'spiderplugin',
        'stlplugin',
        'tinkerplugin',
        'uhbdplugin',
        'xbgfplugin',
        'xsfplugin',
        'xyzplugin' ]

    plugins = [ ]

    clean_re = re.compile("\/\/.*$")
    api_re = re.compile("VMDPLUGIN_API")

    for pref in src_list:
        in_file = glob(molfile_src_path+"/"+pref+".[cC]")[0]
        input = open(in_file).readlines()

        out_file = "src/"+pref+".c"
        plugins.append(pref+".o")
        if in_file[-1:]=='C':
            out_file = out_file + "pp"
            # fix the extern
            input = map(lambda x,c=api_re:c.sub("VMDPLUGIN_EXTERN",x),input)        
        else:
            # get rid of non-ansi C comments
            input = map(lambda x,c=clean_re:c.sub("\n",x),input)
            
        g=open(out_file,'w')
        g.write("/* MACHINE GENERATED FILE, DO NOT EDIT! */\n\n")
        g.write("#define VMDPLUGIN molfile_%s\n"%pref)
        g.write("#define STATIC_PLUGIN 1\n\n")    
        g.write(string.join(input,''))
        g.close()
        g = open("src/objects.make",'w')
        g.write("OBJS="+string.join(plugins,' ')+"\n")
        g.close()


    g=open("src/PlugIOManagerInit.c",'w');
    g.write("/* MACHINE GENERATED FILE, DO NOT EDIT! */\n\n")
    g.write('#include "vmdplugin.h"\n\n')
    g.write('typedef struct _PyMOLGlobals PyMOLGlobals;\n');
    g.write('/* prototypes */')
    for pref in src_list:
        g.write("int molfile_%s_init(void);\n"%pref)
        g.write("int molfile_%s_register(void *,vmdplugin_register_cb *);\n"%pref)
        g.write("int molfile_%s_fini(void);\n"%pref)
    g.write('''

    int PlugIOManagerRegister(PyMOLGlobals *G, vmdplugin_t *);

    int PlugIOManagerInitAll(PyMOLGlobals *G);

    int PlugIOManagerInitAll(PyMOLGlobals *G)
    {
       int ok=1;
    ''')
    for pref in src_list:
        g.write("if(ok) ok = ok && (molfile_%s_init() == VMDPLUGIN_SUCCESS);\n"%pref)
    g.write('''
       if(ok) {
    ''')
    for pref in src_list:
        g.write("if(ok) ok = ok && (molfile_%s_register(G,(vmdplugin_register_cb*)PlugIOManagerRegister) == VMDPLUGIN_SUCCESS);\n"%pref)
    g.write('''
       }
       return ok;
    }
    ''')
    g.write('''

    int PlugIOManagerRegister(PyMOLGlobals *G, vmdplugin_t *);

    int PlugIOManagerFreeAll(void);
    int PlugIOManagerFreeAll(void)
    {
       int ok=1;
    ''')
    for pref in src_list:
        g.write("if(ok) ok = ok && (molfile_%s_fini() == VMDPLUGIN_SUCCESS);\n"%pref)

    g.write('''
       return ok;
    }
    ''')



    os.system("/bin/cp %s/*.h src/"%molfile_src_path)
    os.system("/bin/cp %s/hash.c src/"%molfile_src_path)

