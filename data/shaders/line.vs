#ifndef PYMOL_WEBGL_IOS
#ifdef gl_VertexID_enabled
#extension GL_EXT_gpu_shader4 : require
#endif
#endif

#include webgl_header.vs

attribute vec4 a_Vertex;
attribute vec3 a_Normal;
attribute vec4 a_Color;
attribute float a_interpolate;
#ifndef gl_VertexID_enabled
attribute float a_line_position;
#endif

//#ifdef use_geometry_shaders

#define a_VERTEX gl_Position
varying vec3 a_NORMAL ;
varying vec4 a_COLOR ;
varying vec4 a_COLOR2 ;
varying vec4 a_COLOR_INTERP ;
varying float a_INTERPOLATE ;
varying float a_LINE_POSITION ;
varying float a_fog;
#define bgTextureLookup a_bgTextureLookup
varying vec2 a_bgTextureLookup;

uniform bool isPicking;

void main()
{
#ifdef PYMOL_WEBGL_IOS
  gl_PointSize = g_PointSize;
#endif

  vec3 eye_pos = vec3(g_ModelViewMatrix * a_Vertex);
  a_INTERPOLATE = isPicking ? 0. : a_interpolate;
#ifdef gl_VertexID_enabled
  a_LINE_POSITION = mod(float(gl_VertexID), 2.);
#else
  a_LINE_POSITION = a_line_position;
#endif
  a_VERTEX = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
  a_NORMAL = normalize(g_NormalMatrix * a_Normal);
  a_COLOR_INTERP = a_Color;
  a_COLOR = (1.-a_LINE_POSITION) * a_Color ;
  a_COLOR2 = a_LINE_POSITION * a_Color;
  a_fog = (g_Fog_end + eye_pos.z) * g_Fog_scale;
  a_bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
