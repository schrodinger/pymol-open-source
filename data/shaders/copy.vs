#ifndef PYMOL_WEBGL
#version 120
#endif

uniform sampler2D colorTex;

attribute vec3 a_Vertex;
varying vec2 texcoordAttr;

void main()
{
	vec2 vert = vec2((1. + a_Vertex.xy) / 2.);
	texcoordAttr = vert;
	gl_Position = vec4(a_Vertex.x, a_Vertex.y, 0., 1.);
}
