#include webgl_header.vs

attribute vec3 a_Vertex;
attribute vec3 a_OtherVertex;
attribute float a_UV;
attribute vec4 a_Color;
attribute vec4 a_Color2;
attribute float a_interpolate;

uniform vec2 inv_dimensions;

uniform float line_width;
uniform bool isPicking;

varying vec4 COLOR;
varying vec4 COLOR2;
varying float a_INTERPOLATE;
varying float fog;
varying vec2 bgTextureLookup;

varying float centerdist;
varying float whichEnd;

void main(void) {

  vec2 UV;
  UV.y = mod(a_UV, 2.0);
  UV.x = (a_UV - UV.y) / 2.;
  UV.y = 2. * UV.y - 1.;  // -1. or 1.
  UV.x = 2. * UV.x - 1.;  // -1. or 1.

  float swapPoints = step(0., UV.x);
  vec3 b_Vertex = mix(a_Vertex, a_OtherVertex, swapPoints);
  vec3 b_OtherVertex = mix(a_OtherVertex, a_Vertex, swapPoints);
  whichEnd = swapPoints;

  a_INTERPOLATE = isPicking ? 0. : a_interpolate;
  COLOR = a_Color;
  COLOR2 = a_Color2;

  vec4 eye_pos = g_ModelViewMatrix * vec4(b_Vertex, 1.0);
  vec4 pointA = g_ProjectionMatrix * g_ModelViewMatrix * vec4(b_Vertex, 1.0);
  vec4 pointB = g_ProjectionMatrix * g_ModelViewMatrix * vec4(b_OtherVertex, 1.0);
  pointA.xyz = pointA.xyz / abs(pointA.w);
  pointB.xyz = pointB.xyz / abs(pointB.w);
  pointA.w = 1.;
  pointB.w = 1.;

  // vector perpendicular to centerline
  vec2 perpAB = normalize((pointA.yx - pointB.yx) * inv_dimensions)
      * vec2(1., -1.);

  float width = line_width;

#ifdef line_smooth
  // signed distance from centerline for line smoothing
  centerdist = UV.x * UV.y * line_width;
  if (!isPicking)
    width += 1.0;
  else
    width = max(1.0, width);
#else
  width = max(1.0, width);
#endif

  // offsetting to left or right of centerline
  pointA.xy += width * perpAB * UV.y * inv_dimensions;

  fog = (g_Fog_end + eye_pos.z) * g_Fog_scale;
  gl_Position = pointA;

  bgTextureLookup = (gl_Position.xy / gl_Position.w) / 2.0 + 0.5;
}
