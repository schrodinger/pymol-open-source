#ifndef _CE_TYPES_H
#define _CE_TYPES_H

#include"os_python.h"

/*
// Typical XYZ point and array of points
*/
typedef struct {
	double x;
	double y;
	double z;
} cePoint, *pcePoint;

/*
// An AFP (aligned fragment pair), and list/pointer
*/
typedef struct {
	int first;
	int second;
} afp, *path, **pathCache;

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
pathCache findPath(double** S, double** dA, double**dB, int lenA, int lenB, float D0, float D1, int winSize, int gapMax, int* bufferSize);

// filter through the results and find the best
PyObject* findBest( pcePoint coordsA, pcePoint coordsB, pathCache paths, int bufferSize, int smaller, int winSize );

#endif
