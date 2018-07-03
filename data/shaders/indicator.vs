#ifdef PYMOL_WEBGL_IOS
precision highp float;
#else
#version 120
#endif

#include webgl_header.vs

attribute vec4 a_Vertex;
attribute vec4 a_Color;

varying vec4 COLOR;
varying vec4 POS;

uniform float g_pointSize;

void main()
{
  COLOR = a_Color;
  gl_PointSize = g_pointSize;
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
  POS = gl_Position;
}
