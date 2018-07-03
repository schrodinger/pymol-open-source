#ifndef PYMOL_WEBGL
#version 120
#else
precision highp float;
#endif

uniform sampler2D colorTex;

varying vec2 texcoordAttr;

void main()
{
	gl_FragColor = texture2D(colorTex, texcoordAttr);
}
