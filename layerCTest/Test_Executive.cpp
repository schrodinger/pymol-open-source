#include "Test.h"

#include "Executive.h"
#include "ObjectVolume.h"

using namespace pymol;

#define TEST_SETUP       \
    PyMOLInstance pymol; \
    auto G = pymol.G();  \
    [[maybe_unused]]     \
    bool quiet = true;

#define TEST_SETUP_OBJ_ATOMLESS                  \
    TEST_SETUP                                   \
    auto obj = new ObjectMolecule(G, false);     \
    obj->setName("M1");                          \
    ExecutiveManageObject(G, obj, false, quiet);


#define TEST_SETUP_OBJ                                               \
    TEST_SETUP                                                       \
    ExecutivePseudoatom(G, "M1", "", "PS1", "PSD", "1", "P", "PSDO", \
        "PS", -1.0f, 1, 0.0, 0.0, "", nullptr, -1, -3, 2, 1);        \
    auto obj = ExecutiveFindObjectByName(G, "M1");

// TEST_CASE("ExecutiveManageObject", "[Executive]")
// {
//   TEST_SETUP
//   auto obj = new ObjectMolecule(G, false);
//   obj->setName("M1");
//   ExecutiveManageObject(G, obj, false, quiet);
//   REQUIRE(true);
// }

// TEST_CASE("ExecutiveFindObjectByName", "[Executive]")
// {
//   TEST_SETUP_OBJ_ATOMLESS
//   auto find = ExecutiveFindObjectByName(G, obj->Name);
//   REQUIRE(find == obj);
// }

// TEST_CASE("ExecutiveFindObject", "[Executive]")
// {
//   TEST_SETUP_OBJ_ATOMLESS
//   auto find = ExecutiveFindObject<ObjectMolecule>(G, obj->Name);
//   REQUIRE(find == obj);
//   auto bad_find = ExecutiveFindObject<ObjectVolume>(G, obj->Name);
//   REQUIRE(bad_find == nullptr);
// }

// TEST_CASE("ExecutiveGetNames", "[Executive]")
// {
//   TEST_SETUP_OBJ
//   auto type = 1; //objects
//   auto enabled_only = false;
//   auto selection = obj->Name;
//   auto names_result = ExecutiveGetNames(G, type, enabled_only, selection);
//   const auto& names = *names_result;
//   REQUIRE(names.size() == 1u);
//   REQUIRE(pymol::zstring_view(names[0]) == obj->Name);
// }
