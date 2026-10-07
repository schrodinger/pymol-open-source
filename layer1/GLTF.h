/*
 * PyMOL glTF 2.0 / GLB export
 *
 * (c) 2026 Schrodinger, Inc.
 */

#pragma once

#include "Base.h"

#ifdef _HAVE_JSON

/**
 * Render the current scene as a GLB (binary glTF 2.0) file. Expands all
 * objects into ray-tracing primitives, tessellates spheres/cylinders/cones
 * into triangle meshes, groups by transparency, and writes the GLB container
 * (header + JSON chunk + BIN chunk) into the output VLA.
 * @param I ray context with primitive data
 * @param[in,out] vla_ptr pre-allocated VLA; on return, resized to hold the
 *   GLB binary data (little endian), or resized to 0 if the scene has no
 *   exportable geometry or does not fit the 32-bit GLB container limit
 * @note The scene is oriented like the current view, unless
 *   geometry_export_mode is 1, which keeps model space coordinates.
 * @note Text label and ellipsoid primitives are not supported yet, they are
 *   skipped with a warning.
 */
void RayRenderGLB(CRay* I, char** vla_ptr);

#endif
