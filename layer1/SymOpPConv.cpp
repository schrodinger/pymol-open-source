/**
 * @file Python serialization of SymOp
 *
 * (c) Schrodinger, Inc.
 */

#include "SymOpPConv.h"
#include "SymOp.h"

PyObject* PConvToPyObject(pymol::SymOp const& symop)
{
  std::string buffer;
  if (symop) {
    buffer = symop.to_string();
  }
  return PyString_FromString(buffer.c_str());
}

bool PConvFromPyObject(PyMOLGlobals*, PyObject* obj, pymol::SymOp& out)
{
  auto str = PyString_AsSomeString(obj);
  out.reset(str.c_str());
  return true;
}
