#include webgl_header.vs

attribute vec4 a_Vertex;
attribute vec3 a_Normal;
attribute vec4 a_Color;

varying vec3 packed_data_0 ;

//varying vec3 N;
#define NORMAL  packed_data_0

#ifdef PURE_OPENGL_ES_2
varying vec4 COLOR;
#else
// using the built-in allows to use glShadeModel
#define COLOR gl_FrontColor = gl_BackColor
#endif

varying float fog;
varying vec2 bgTextureLookup;

void main()
{
#ifdef PYMOL_WEBGL_IOS
  gl_PointSize = g_PointSize;
#endif

  NORMAL = normalize(g_NormalMatrix * a_Normal);
  vec3 eye_pos = vec3(g_ModelViewMatrix * a_Vertex);

  COLOR = a_Color;

  fog = (g_Fog_end + eye_pos.z) * g_Fog_scale;
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
