#include <stdio.h>

#ifndef _WIN32
#include <stdlib.h>
#include <unistd.h>
#endif

#include <fstream>
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

namespace pymol {
namespace test {

TmpFILE::TmpFILE()
{
#ifdef _WIN32
    tmpFilename.resize(L_tmpnam_s);
    tmpnam_s(&tmpFilename[0], tmpFilename.size());
    tmpFilename.resize(strlen(tmpFilename.c_str()));

    // file 'touch'
    std::ofstream(tmpFilename);
#else
    tmpFilename = P_tmpdir;

    if (!tmpFilename.empty() && tmpFilename.back() != '/') {
      tmpFilename += '/';
    }

    tmpFilename.append("tmppymoltestXXXXXX");

    close(mkstemp(&tmpFilename[0]));
#endif
}

} // namespace test
} // namespace pymol
