# get the function we need from deep inside PyMOL

from pymol.wizard.annotation import load_annotated_sdf

# add this function to the command language

cmd.extend('load_annotated_sdf',load_annotated_sdf)

# activate the annotation wizard

wizard annotation

# make the prompt transparent 

set wizard_prompt_mode,3

# now load an SD file

load_annotated_sdf demo.sdf

