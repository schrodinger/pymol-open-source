#include webgl_header.fs

#ifdef PYMOL_WEBGL_IOS
precision highp float;
#else
#version 120
#endif

uniform sampler2D textureMap;
uniform vec2 textureLookup ;
uniform vec2 textureScale ;
uniform float g_pointSize;
uniform vec4 viewport;

#include anaglyph_header.fs

varying vec4 POS;

/*
* On graphics drivers which are not buggy, this should be
* equivalent to gl_PointCoord.
*/
vec2 get_pointcoord() {
 vec2 p = viewport.xy + 0.5 * (POS.xy / POS.w + 1.0) * viewport.zw;
 return (gl_FragCoord.xy - p) / g_pointSize + 0.5;
}

void main()
{
 gl_FragColor = texture2D(textureMap, textureLookup + get_pointcoord() * textureScale);
 PostLightingEffects();
}
