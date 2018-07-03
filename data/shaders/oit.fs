#ifndef PYMOL_WEBGL
#version 120
#else
precision highp float;
#endif

uniform sampler2D accumTex;
uniform sampler2D revealageTex;
uniform float isRight;

varying vec2 texcoordAttr;

#define texture texture2D

void main()
{
        vec4 accum = texture(accumTex, texcoordAttr).rgba;
        float r = accum.a;
#ifdef ANAGLYPH
        // need to set two components because anaglyph masks out colors for each eye
#ifdef GLEW_VERSION_3_0
        accum.a = texture(revealageTex, texcoordAttr).r;
#else
        vec2 sum = texture(revealageTex, texcoordAttr).rg;
        accum.a = (1.-isRight) * sum.x + isRight * sum.y;
#endif

#else
        accum.a = texture(revealageTex, texcoordAttr).r;
#endif
	gl_FragColor = vec4(accum.rgb / clamp(accum.a, 1e-4, 5e4), r);
}
