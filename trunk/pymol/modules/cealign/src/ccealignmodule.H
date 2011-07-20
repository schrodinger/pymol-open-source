//
// ccealign -- structural alignment plugin for PyMol
//
//////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) Jason Vertrees.
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//      * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//      * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in
//      the documentation and/or other materials provided with the
//      distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
//  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
//  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//////////////////////////////////////////////////////////////////////////////
#include"os_python.h"

// include the template numerical toolkit
#include "tnt/tnt.h"
// for reflections
#include "jama_lu.h"
// for the SVD
#include "jama_svd.h"

/////////////////////////////////////////////////////////////////////////////
// Structs, typedefs, declarations that need C++ linkage
/////////////////////////////////////////////////////////////////////////////
#define TA1 TNT::Array1D
#define TA2 TNT::Array2D
// tranpose of a matrix
TA2<double> transpose( const TA2<double>& v );

// C-linkage declarations
#ifdef __cplusplus
extern "C" { 
#endif

#include"ce_types.h"
#include"os_std.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////////
// Function Declarations
/////////////////////////////////////////////////////////////////////////////
// Calculates the CE Similarity Matrix
double** calcS(double** d1, double** d2, int lenA, int lenB, int wSize);

// calculates a simple distance matrix
double** calcDM(pcePoint coords, int len);

// Converter: Python Object -> C Structs
pcePoint getCoords( PyObject* L, int len );

// Optimal path finding algorithm (CE).
pathCache findPath(double** S, double** dA, double**dB, int lenA, int lenB, int winSize, int* bufferSize);

// filter through the results and find the best
PyObject* findBest( pcePoint coordsA, pcePoint coordsB, pathCache paths, int bufferSize, int smaller, int winSize );

#ifdef __cplusplus
}
#endif
