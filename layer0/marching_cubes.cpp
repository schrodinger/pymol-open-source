/**
 * @file Marching cubes
 *
 * Based on https://github.com/ilastik/marching_cubes
 *
 * License: BSD-3-Clause
 *
 * Copyright 2018, The ilastik development team
 * Copyright 2020, Schrodinger, Inc.
 */

#include "marching_cubes.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <vector>

#ifdef PYMOL_OPENMP
#include <omp.h>
#endif

namespace mc
{

static const unsigned short EDGE_TABLE[256] = {0x0, 0x109, 0x203, 0x30a, 0x406,
    0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f,
    0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f,
    0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0,
    0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6,
    0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265,
    0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9,
    0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa,
    0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3,
    0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3,
    0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,
    0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0, 0x950, 0x859, 0xb53, 0xa5a,
    0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759,
    0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5,
    0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0xf66,
    0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af,
    0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0, 0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f,
    0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230, 0xe90,
    0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596,
    0x29a, 0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905,
    0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0};

static const signed char TRIANGLE_TABLE[256][16] = {
    {-1},
    {0, 8, 3, -1},
    {0, 1, 9, -1},
    {1, 8, 3, 9, 8, 1, -1},
    {1, 2, 10, -1},
    {0, 8, 3, 1, 2, 10, -1},
    {9, 2, 10, 0, 2, 9, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1},
    {3, 11, 2, -1},
    {0, 11, 2, 8, 11, 0, -1},
    {1, 9, 0, 2, 3, 11, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1},
    {3, 10, 1, 11, 10, 3, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1},
    {9, 8, 10, 10, 8, 11, -1},
    {4, 7, 8, -1},
    {4, 3, 0, 7, 3, 4, -1},
    {0, 1, 9, 8, 4, 7, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1},
    {1, 2, 10, 8, 4, 7, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1},
    {8, 4, 7, 3, 11, 2, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1},
    {9, 5, 4, -1},
    {9, 5, 4, 0, 8, 3, -1},
    {0, 5, 4, 1, 5, 0, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1},
    {1, 2, 10, 9, 5, 4, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1},
    {9, 5, 4, 2, 3, 11, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1},
    {9, 7, 8, 5, 7, 9, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1},
    {1, 5, 3, 3, 5, 7, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1},
    {10, 6, 5, -1},
    {0, 8, 3, 5, 10, 6, -1},
    {9, 0, 1, 5, 10, 6, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1},
    {1, 6, 5, 2, 6, 1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1},
    {2, 3, 11, 10, 6, 5, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1},
    {5, 10, 6, 4, 7, 8, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1},
    {10, 4, 9, 6, 4, 10, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1},
    {0, 2, 4, 4, 2, 6, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1},
    {6, 4, 8, 11, 6, 8, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1},
    {7, 3, 2, 6, 7, 2, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1},
    {7, 11, 6, -1},
    {7, 6, 11, -1},
    {3, 0, 8, 11, 7, 6, -1},
    {0, 1, 9, 11, 7, 6, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1},
    {10, 1, 2, 6, 11, 7, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1},
    {7, 2, 3, 6, 2, 7, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1},
    {6, 8, 4, 11, 8, 6, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1},
    {0, 4, 2, 4, 6, 2, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1},
    {4, 9, 5, 7, 6, 11, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1},
    {0, 3, 8, 5, 6, 10, -1},
    {10, 5, 6, -1},
    {11, 5, 10, 7, 5, 11, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1},
    {9, 8, 7, 5, 9, 7, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1},
    {0, 4, 5, 1, 0, 5, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1},
    {9, 4, 5, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1},
    {4, 0, 3, 7, 4, 3, -1},
    {4, 8, 7, -1},
    {9, 10, 8, 10, 11, 8, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1},
    {3, 1, 10, 11, 3, 10, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1},
    {0, 2, 11, 8, 0, 11, -1},
    {3, 2, 11, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1},
    {9, 10, 2, 0, 9, 2, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1},
    {1, 10, 2, -1},
    {1, 3, 8, 9, 1, 8, -1},
    {0, 9, 1, -1},
    {0, 3, 8, -1},
    {-1},
};

constexpr unsigned EDGES_PER_CELL = 3;

static size_t vertexId(size_t x, size_t y, size_t z, size_t xDim, size_t yDim)
{
  return EDGES_PER_CELL * (z * yDim * xDim + y * xDim + x);
}

inline size_t edgeId2z(size_t eid, size_t xDim, size_t yDim)
{
  return (eid / EDGES_PER_CELL) / (yDim * xDim);
}

static size_t edgeId(
    size_t x, size_t y, size_t z, size_t edgeNumber, size_t xDim, size_t yDim)
{
  switch (edgeNumber) {
  case 0:
    return vertexId(x, y, z, xDim, yDim) + 1;
  case 1:
    return vertexId(x, y + 1, z, xDim, yDim);
  case 2:
    return vertexId(x + 1, y, z, xDim, yDim) + 1;
  case 3:
    return vertexId(x, y, z, xDim, yDim);
  case 4:
    return vertexId(x, y, z + 1, xDim, yDim) + 1;
  case 5:
    return vertexId(x, y + 1, z + 1, xDim, yDim);
  case 6:
    return vertexId(x + 1, y, z + 1, xDim, yDim) + 1;
  case 7:
    return vertexId(x, y, z + 1, xDim, yDim);
  case 8:
    return vertexId(x, y, z, xDim, yDim) + 2;
  case 9:
    return vertexId(x, y + 1, z, xDim, yDim) + 2;
  case 10:
    return vertexId(x + 1, y + 1, z, xDim, yDim) + 2;
  case 11:
    return vertexId(x + 1, y, z, xDim, yDim) + 2;
  default:
    assert(false);
    return -1;
  }
}

inline float mix(float v0, float v1, float a)
{
  return v0 * (1 - a) + v1 * a;
}

inline Point mix(const Point& v0, const Point& v1, float a)
{
  return {
      mix(v0.x, v1.x, a),
      mix(v0.y, v1.y, a),
      mix(v0.z, v1.z, a),
  };
}

static Point normalize(Point p)
{
  auto length = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  if (length > 0.) {
    p[0] /= length;
    p[1] /= length;
    p[2] /= length;
  } else {
    p[0] = 1;
  }
  return p;
}

/**
 * Calculate the position and (optionally) the normal on the given edge at the
 * given isolevel.
 */
static Point calculateIntersection(const Field& volume, float isoLevel,
    size_t x, size_t y, size_t z, size_t edgeNumber, Point* normal = nullptr)
{
  auto x1 = x;
  auto x2 = x;
  auto y1 = y;
  auto y2 = y;
  auto z1 = z;
  auto z2 = z;

  switch (edgeNumber) {
  case 0:
    y2 += 1;
    break;
  case 1:
    y1 += 1;
    x2 += 1;
    y2 += 1;
    break;
  case 2:
    x1 += 1;
    y1 += 1;
    x2 += 1;
    break;
  case 3:
    x1 += 1;
    break;
  case 4:
    z1 += 1;
    y2 += 1;
    z2 += 1;
    break;
  case 5:
    y1 += 1;
    z1 += 1;
    x2 += 1;
    y2 += 1;
    z2 += 1;
    break;
  case 6:
    x1 += 1;
    y1 += 1;
    z1 += 1;
    x2 += 1;
    z2 += 1;
    break;
  case 7:
    x1 += 1;
    z1 += 1;
    z2 += 1;
    break;
  case 8:
    z2 += 1;
    break;
  case 9:
    y1 += 1;
    y2 += 1;
    z2 += 1;
    break;
  case 10:
    x1 += 1;
    y1 += 1;
    x2 += 1;
    y2 += 1;
    z2 += 1;
    break;
  case 11:
    x1 += 1;
    x2 += 1;
    z2 += 1;
    break;
  default:
    assert(false);
  }

  auto const pos1 = volume.get_point(x1, y1, z1);
  auto const pos2 = volume.get_point(x2, y2, z2);
  auto const val1 = volume.get(x1, y1, z1);
  auto const val2 = volume.get(x2, y2, z2);
  assert(val1 != val2);
  auto const frac = (isoLevel - val1) / (val2 - val1);

  if (normal) {
    auto grad1 = volume.get_gradient(x1, y1, z1);
    auto grad2 = volume.get_gradient(x2, y2, z2);
    *normal = normalize(mix(grad1, grad2, frac));
  }

  return mix(pos1, pos2, frac);
}

struct IdPoint {
  size_t number = 0;
  Point point;
  Point normal;
};

struct Triangle {
  size_t pointId[3];
};

/**
 * The marching cubes algorithm as described here:
 * http://paulbourke.net/geometry/polygonise/
 *
 * @param volume The data field
 * @param isoLevel Minimum isoLevel, all values >= isoLevel will contribute to
 * the mesh
 * @param gradient_normals Compute normals based on field gradient. If false,
 * then don't compute normals.
 * @return The iso-surface mesh
 */
Mesh march(const Field& volume, float isoLevel, bool gradient_normals)
{
  auto const xDim = volume.xDim();
  auto const yDim = volume.yDim();
  auto const zDim = volume.zDim();

#ifndef PYMOL_OPENMP
  using isocheck_bool_t = bool;
#else
  // writes to std::vector<bool> are not thread-safe
  using isocheck_bool_t = char;
#endif

  // pre-compute isovalue check for better performance
  std::vector<isocheck_bool_t> isocheck(xDim * yDim * zDim);

#pragma omp parallel for
  for (int z = 0; z < zDim; ++z) {
    for (size_t y = 0; y < yDim; ++y) {
      auto const offset = xDim * y + xDim * yDim * z;
      for (size_t x = 0; x < xDim; ++x) {
        isocheck[x + offset] = volume.get(x, y, z) < isoLevel;
      }
    }
  }

  auto const get_isocheck = [&](size_t x, size_t y, size_t z) -> bool {
    return isocheck[x + xDim * y + xDim * yDim * z];
  };

  auto const xEnd = xDim - 1;
  auto const yEnd = yDim - 1;
  auto const zEnd = zDim - 1;

  // With OpenMP, we use one triangles vector per thread, and one vertexMap per
  // z-index, and let OpenMP distribute the z-index runs across threads.
  std::vector<std::vector<Triangle>> trianglesVec(1);
  std::vector<std::unordered_map<size_t, IdPoint>> vertexMapVec(1);

#ifndef PYMOL_OPENMP
#define omp_get_thread_num() 0
#define vertexMappingGet(eid) vertexMapVec[0][eid]
#else
#define vertexMappingGet(eid) vertexMapVec[edgeId2z(eid, xDim, yDim)][eid]

  trianglesVec.resize(omp_get_max_threads());
  vertexMapVec.resize(zDim);

#pragma omp parallel for
#endif
  for (int z = 0; z < zEnd; ++z) {
    auto& triangles = trianglesVec[omp_get_thread_num()];

    for (size_t y = 0; y < yEnd; ++y) {
      for (size_t x = 0; x < xEnd; ++x) {
        size_t tableIndex = 0;
        if (get_isocheck(x, y, z))
          tableIndex |= 1;
        if (get_isocheck(x, y + 1, z))
          tableIndex |= 2;
        if (get_isocheck(x + 1, y + 1, z))
          tableIndex |= 4;
        if (get_isocheck(x + 1, y, z))
          tableIndex |= 8;
        if (get_isocheck(x, y, z + 1))
          tableIndex |= 16;
        if (get_isocheck(x, y + 1, z + 1))
          tableIndex |= 32;
        if (get_isocheck(x + 1, y + 1, z + 1))
          tableIndex |= 64;
        if (get_isocheck(x + 1, y, z + 1))
          tableIndex |= 128;

        if (EDGE_TABLE[tableIndex] == 0) {
          continue;
        }

        auto const storeEdgeVertex = [&](size_t edgeNumber) {
          if ((EDGE_TABLE[tableIndex] & (1 << edgeNumber)) != 0) {
            auto eid = edgeId(x, y, z, edgeNumber, xDim, yDim);
            auto& edge = vertexMappingGet(eid);
            edge.point = calculateIntersection(volume, isoLevel, x, y, z,
                edgeNumber, gradient_normals ? &edge.normal : nullptr);
          }
        };

        storeEdgeVertex(3);
        storeEdgeVertex(0);
        storeEdgeVertex(8);

        if (x == xEnd - 1) {
          storeEdgeVertex(2);
          storeEdgeVertex(11);

          if (y == yEnd - 1) {
            storeEdgeVertex(10);
          }
        }

        if (y == yEnd - 1) {
          storeEdgeVertex(1);
          storeEdgeVertex(9);

          if (z == zEnd - 1) {
            storeEdgeVertex(5);
          }
        }

        if (z == zEnd - 1) {
          storeEdgeVertex(4);
          storeEdgeVertex(7);

          if (x == xEnd - 1) {
            storeEdgeVertex(6);
          }
        }

        auto const* tri_table_row = TRIANGLE_TABLE[tableIndex];
        for (size_t i = 0; tri_table_row[i] != -1; i += 3) {
          auto pointId0 = edgeId(x, y, z, tri_table_row[i], xDim, yDim);
          auto pointId1 = edgeId(x, y, z, tri_table_row[i + 1], xDim, yDim);
          auto pointId2 = edgeId(x, y, z, tri_table_row[i + 2], xDim, yDim);
          triangles.push_back({pointId0, pointId1, pointId2});
        }
      }
    }
  }

  Mesh mesh;
  for (auto const& vertexMap : vertexMapVec) {
    mesh.vertexCount += vertexMap.size();
  }
  for (auto const& triangles : trianglesVec) {
    mesh.faceCount += triangles.size();
  }

  mesh.faces.reset(new size_t[mesh.faceCount * 3]);
  mesh.vertices.reset(new Point[mesh.vertexCount]);
  if (gradient_normals) {
    mesh.normals.reset(new Point[mesh.vertexCount]);
  }

  size_t index = 0;
  for (auto& vertexMap : vertexMapVec) {
    for (auto& pair : vertexMap) {
      mesh.vertices[index] = pair.second.point;
      if (gradient_normals) {
        mesh.normals[index] = pair.second.normal;
      }

      // index for faces
      pair.second.number = index++;
    }
  }

  index = 0;
  for (auto const& triangles : trianglesVec) {
    for (const auto& triangle : triangles) {
      mesh.faces[index++] = vertexMappingGet(triangle.pointId[0]).number;
      mesh.faces[index++] = vertexMappingGet(triangle.pointId[1]).number;
      mesh.faces[index++] = vertexMappingGet(triangle.pointId[2]).number;
    }
  }

  return mesh;
}

/**
 * Calculate triangle-based normals
 */
void calculateNormals(Mesh& mesh)
{
  size_t vertexCount = mesh.vertexCount;
  const Point* vertices = mesh.vertices.get();
  size_t triangleCount = mesh.faceCount;
  const size_t* triangles = mesh.faces.get();

  mesh.normals.reset(new Point[vertexCount]);
  auto normals = mesh.normals.get();

#pragma omp parallel for
  for (int i = 0; i < vertexCount; ++i) {
    normals[i] = {0, 0, 0};
  }

#pragma omp parallel for
  for (int i = 0; i < triangleCount; ++i) {
    size_t const id0 = triangles[i * 3];
    size_t const id1 = triangles[i * 3 + 1];
    size_t const id2 = triangles[i * 3 + 2];
    Point const vec1{
        vertices[id1][0] - vertices[id0][0],
        vertices[id1][1] - vertices[id0][1],
        vertices[id1][2] - vertices[id0][2],
    };
    Point const vec2{
        vertices[id2][0] - vertices[id0][0],
        vertices[id2][1] - vertices[id0][1],
        vertices[id2][2] - vertices[id0][2],
    };
    Point const normal{
        vec1[2] * vec2[1] - vec1[1] * vec2[2],
        vec1[0] * vec2[2] - vec1[2] * vec2[0],
        vec1[1] * vec2[0] - vec1[0] * vec2[1],
    };
#pragma omp critical
    {
      normals[id0][0] += normal[0];
      normals[id0][1] += normal[1];
      normals[id0][2] += normal[2];
      normals[id1][0] += normal[0];
      normals[id1][1] += normal[1];
      normals[id1][2] += normal[2];
      normals[id2][0] += normal[0];
      normals[id2][1] += normal[1];
      normals[id2][2] += normal[2];
    }
  }

#pragma omp parallel for
  for (int i = 0; i < vertexCount; ++i) {
    normals[i] = normalize(normals[i]);
  }
}

Point Field::get_gradient(size_t x, size_t y, size_t z) const
{
  size_t xx[] = {x == 0 ? x : x - 1, std::min(x + 1, xDim() - 1)};
  size_t yy[] = {y == 0 ? y : y - 1, std::min(y + 1, yDim() - 1)};
  size_t zz[] = {z == 0 ? z : z - 1, std::min(z + 1, zDim() - 1)};
  return {
      (get(xx[0], y, z) - get(xx[1], y, z)) / std::max<int>(1, xx[1] - xx[0]),
      (get(x, yy[0], z) - get(x, yy[1], z)) / std::max<int>(1, yy[1] - yy[0]),
      (get(x, y, zz[0]) - get(x, y, zz[1])) / std::max<int>(1, zz[1] - zz[0]),
  };
}

} // namespace mc
