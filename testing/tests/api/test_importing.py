from pymol import cmd
from pymol import test_utils
from pymol.querying import cif_get_array


@test_utils.requires_version("3.0")
def test_bcif():
    cmd.load(test_utils.datafile("115d.bcif.gz"))
    assert cmd.count_atoms() == 407

@test_utils.requires_version("3.0")
def test_bcif_array():
    obj_name = "foo"
    cmd.set('cif_keepinmemory', 1)
    cmd.load(test_utils.datafile("115d.bcif.gz"), object=obj_name)
    arr = cif_get_array(obj_name, "_pdbx_database_status.entry_id", "s")
    assert arr == ["115D"]

    arr = cif_get_array(obj_name, "_entity_poly.pdbx_strand_id", "s")
    assert arr == ["A,B"]

    arr = cif_get_array(obj_name, "_pdbx_struct_oper_list.name", "s")
    assert arr == ["1_555"]

    arr = cif_get_array(obj_name, "_pdbx_struct_assembly.oligomeric_count", "i")
    assert arr == [2]
