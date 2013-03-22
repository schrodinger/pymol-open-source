#ifdef PYMOL_IOS
precision highp float;
#else
#version 120
#endif

attribute vec4 a_Vertex;
attribute vec4 a_Color;

varying vec4 COLOR;

#ifdef PYMOL_IOS
uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
#endif

uniform float g_pointSize;

void main()
{
  COLOR = a_Color;
  gl_PointSize = g_pointSize;
#ifdef PYMOL_IOS
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
#else
  gl_Position = gl_ModelViewProjectionMatrix * a_Vertex;
#endif
}
