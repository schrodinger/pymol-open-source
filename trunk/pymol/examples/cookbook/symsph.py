# symshp: create a symmetry-expanded sphere about a selection

# usage:
#
#    symexp name [,selection [,cutoff ]]

from pymol import cmd as global_cmd

def symsph(name, selection="sele", cutoff=20.0, self_cmd=global_cmd):
    cutoff = float(cutoff)
    prefix = selection+"_symarea_"
    tmp_obj = selection+"_tmp"
    if selection not in self_cmd.get_names("selections"):
        print " error: '"+selection+"' is not defined."
        return self_cmd.DEFAULT_FAILURE
    if not self_cmd.count_atoms(selection):
        print " error: '"+selection+"' contains no atoms."
        return self_cmd.DEFAULT_FAILURE
    obj_list = self_cmd.get_object_list(selection)
    if len(obj_list)!=1:
        print script_name+" error: '"+selection+"' must only span one object.'"
        return self_cmd.DEFAULT_FAILURE
    obj = obj_list[0]
    cmd.center(selection)
    cmd.pseudoatom(tmp_obj)
    cmd.delete(prefix+"*")
    cmd.symexp(prefix,obj,tmp_obj,cutoff,segi=1)
    cmd.create(name,"("+obj+" or "+prefix+"*) within %1.9f of %s"%(cutoff,tmp_obj))
    cmd.delete(tmp_obj)
    cmd.delete(prefix+"*")
    
#symsph("sele",20)

cmd.extend("symsph",symsph)


    

