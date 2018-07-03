#include webgl_header.vs

// Sphere impostor vertex shader

uniform float sphere_size_scale;

attribute vec4 a_vertex_radius;
attribute vec4 a_Color;
attribute float a_rightUpFlags;

varying vec4 COLOR;
varying vec3 sphere_center;
varying float radius2;
varying vec3 point;
varying vec2 bgTextureLookup;

/*
 * horizontial and vertical adjustment of outer tangent hitting the
 * impostor quad, in model view space.
 */
vec2 outer_tangent_adjustment(vec3 center, float radius_sq) {
    vec2 xy_dist = vec2(length(center.xz), length(center.yz));

    // without clamping, this caused flickering (divide-by-zero)
    vec2 cos_a = clamp(center.z / xy_dist, -1., 1.);
    vec2 cos_b = xy_dist / sqrt(radius_sq + (xy_dist * xy_dist));

    // numerically more stable version of:
    // vec2 tan_ab_sq = pow(tan(acos(cos_b) + acos(cos_a)), 2);
    vec2 cos_ab = (cos_a * cos_b + sqrt(
                (1. - cos_a * cos_a) *
                (1. - cos_b * cos_b)));
    vec2 cos_ab_sq = cos_ab * cos_ab;
    vec2 tan_ab_sq = (1. - cos_ab_sq) / cos_ab_sq;

    vec2 adjustment = sqrt(tan_ab_sq + 1.);

    // max out (empirical) to avoid exploding (can happen for spheres outside of the viewport)
    return min(adjustment, 10.);
}

void main(void)
{
    // Get billboard attributes
    float radius = a_vertex_radius.w * sphere_size_scale;

    // support uniform scaling
    radius /= length(g_NormalMatrix[0]);

    float right = -1. + 2.*mod(a_rightUpFlags, 2.);
    float up = -1. + 2.*floor(mod(a_rightUpFlags/2., 2.));
    vec4 tmppos = g_ModelViewMatrix * vec4(a_vertex_radius.xyz, 1.);

    COLOR = a_Color;
    radius2 = radius * radius; // compute squared radius 

    // uni-radius corner offset
    vec2 corner_offset = vec2(right, up);

#ifndef ortho
    // horizontial and vertical adjustment due to projection
    corner_offset *= outer_tangent_adjustment(tmppos.xyz, radius2);
#endif

    // corner vertex
    vec4 eye_space_pos = tmppos;
    eye_space_pos.xy += radius * corner_offset;

    // Compute sphere position in modelview space
    sphere_center = tmppos.xyz / tmppos.w;

    // Compute ray direction and origin point
    point = eye_space_pos.xyz / eye_space_pos.w;

#ifndef PYMOL_WEBGL_IOS
    // Pass the transformed vertex for clipping plane calculations
    gl_ClipVertex = eye_space_pos;
#endif

    // Pass transformed vertex
    gl_Position = g_ProjectionMatrix * eye_space_pos;
    bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}

