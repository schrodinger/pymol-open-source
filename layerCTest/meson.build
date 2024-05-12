src_files += files(
    'Test.cpp',
    'Test_Algorithm.cpp',
    'Test_Bezier.cpp',
    'Test_cache_ptr.cpp',
    'Test_CCrystal.cpp',
    'Test_CifFile.cpp',
    'Test_Classic_VLA.cpp',
    'Test_copyable_ptr.cpp',
    'Test_Event.cpp',
    'Test_Executive.cpp',
    'Test_Image.cpp',
    'Test_List.cpp',
    'Test_Picking.cpp',
    'Test_Result.cpp',
    'Test_ScrollBar.cpp',
    'Test_ShaderPreprocessor.cpp',
    'Test_SymOp.cpp',
    'Test_Test.cpp',
    'Test_TTT.cpp',
    'Test_type_traits.cpp',
    'Test_Util.cpp',
    'Test_VLA.cpp',
    'Test_zstring_view.cpp',
)


dependencies += subproject('catch2').get_variable('catch2_dep')
include_dirs += include_directories('.')
macros += '-D_PYMOL_CTEST'

