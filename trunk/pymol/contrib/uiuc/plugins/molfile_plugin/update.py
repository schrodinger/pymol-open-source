if 1:
    print "DANGER DO NOT RUN UNTIL changes to the following files have been backported..."
    print "dcdplugin.c, gridplugin.c, endianswap.h, "
    print "dtrplugin.cpp, fs4plugin.cpp, maeffplugin.cpp"
else:

    import re
    import string
    import os
    from glob import glob

#    molfile_src_path = "/home/warren/06duo/software/vmd/plugins/molfile_plugin/src"
    molfile_src_path = "/Users/delwarl/tmp/plugins/molfile_plugin/src"    

    src_list=[
        'avsplugin',
#        'babelplugin', # requires openbabel
        'basissetplugin', #
        'basissetplugin', #
        'bgfplugin',
        'binposplugin',
        'biomoccaplugin',
        'brixplugin',
        'carplugin',
        'ccp4plugin',
#        'cdfplugin', # requires netcdf
        'corplugin',
        'cpmdlogplugin', #
        'cpmdplugin',
        'crdplugin',
        'cubeplugin',
        'dcdplugin',
        'dlpolyplugin',
        'dsn6plugin',
        'dtrplugin', #
        'dxplugin',
        'edmplugin',
        'fs4plugin',
        'gamessplugin',
        'gaussianplugin', #
        'graspplugin',
        'grdplugin',
        'gridplugin',
        'gromacsplugin',
#        'hoomdplugin', # requires expat
        'jsplugin', #
#        'lammpsplugin', # requires gz
        'maeffplugin', #
        'mapplugin',
        'mdfplugin',
        'mmcif', #
        'mol2plugin',
        'moldenplugin',
        'mrcplugin', #
        'msmsplugin',
        'namdbinplugin',
#        'netcdfplugin', # requires netcdf
        'parm7plugin',
        'parmplugin',
        'pbeqplugin', #
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
        'vaspchgcarplugin', #
        'vaspoutcarplugin', #
        'vaspposcarplugin', #
        'vaspxdatcarplugin', #
        'vaspxmlplugin', #
        'vtfplugin', #
#        'webpdbplugin', # tcl dependent
        'xbgfplugin',
        'xsfplugin',
        'xyzplugin']

    plugins = [ ]

    clean_re = re.compile("\/\/.*$")
    api_re = re.compile("VMDPLUGIN_API")

    for pref in src_list:
        print pref
        in_file = glob(molfile_src_path+"/"+pref+".[cC]*")[0]
        input = open(in_file).readlines()

        out_file = "src/"+pref+".c"
        plugins.append(pref+".o")
        if (in_file[-1:]=='C') or (in_file[-3:]=='cxx'):
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
        output = string.join(input,'')
        output.replace('(vmdplugin_t *)','(vmdplugin_t *)(void*)')
        g.write(output)
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
        g.write("int molfile_%s_register(void *,vmdplugin_register_cb);\n"%pref)
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
        g.write("if(ok) ok = ok && (molfile_%s_register(G,(vmdplugin_register_cb)PlugIOManagerRegister) == VMDPLUGIN_SUCCESS);\n"%pref)
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
    }\n''')


    os.system("/bin/cp %s/*.h src/"%molfile_src_path)
    os.system("/bin/cp %s/hash.c src/"%molfile_src_path)

    os.system("/bin/chmod -x src/*")
