# symshp: create a symmetry-expanded sphere about a selection

from pymol import cmd as global_cmd

def symsph(object_name, target_sele="sele", radius=20.0, self_cmd=global_cmd):
    radius = float(radius)
    prefix = target_sele+"_symarea_"
    tmp_obj = target_sele+"_tmp"
    if target_sele not in self_cmd.get_names("selections"):
        print " error: '"+target_sele+"' is not defined."
        return self_cmd.DEFAULT_FAILURE
    if not self_cmd.count_atoms(target_sele):
        print " error: '"+target_sele+"' contains no atoms."
        return self_cmd.DEFAULT_FAILURE
    obj_list = self_cmd.get_object_list(target_sele)
    if len(obj_list)!=1:
        print script_name+" error: '"+target_sele+"' must only span one object.'"
        return self_cmd.DEFAULT_FAILURE
    obj = obj_list[0]
    cmd.center(target_sele)
    cmd.pseudoatom(tmp_obj)
    cmd.delete(prefix+"*")
    cmd.symexp(prefix,obj,tmp_obj,radius,segi=1)
    cmd.create("symsph","("+obj+" or "+prefix+"*) within %1.9f of %s"%(radius,tmp_obj))
    cmd.delete(tmp_obj)
    cmd.delete(prefix+"*")
    
#symsph("sele",20)

cmd.extend("symsph",symsph)


    

