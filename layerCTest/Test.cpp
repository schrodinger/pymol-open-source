#include <iostream>
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include "Test.h"
#include "TestCmdTest2.h"

using PyMOL_TestAPI = pymol::test::PYMOL_TEST_API;

PyObject *PyMOL_TestAPI::PYMOL_TEST_SUCCESS = PConvAutoNone(Py_None);
PyObject *PyMOL_TestAPI::PYMOL_TEST_FAILURE = Py_BuildValue("i", -1);

PyObject *CmdTest2(PyObject *, PyObject *) {
  int argc = 1;
  char argv0[] = "pymol";
  char *argv[] = {argv0};
  auto result = Catch::Session().run(argc, argv);
  if (!result) {
    return PyMOL_TestAPI::PYMOL_TEST_SUCCESS;
  } else {
    return PyMOL_TestAPI::PYMOL_TEST_FAILURE;
  }
}
