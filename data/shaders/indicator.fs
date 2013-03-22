#ifdef PYMOL_IOS
precision mediump float;
#else
#version 120
#endif

uniform sampler2D textureMap;
uniform vec2 textureLookup ;
uniform vec2 textureScale ;

#include ANAGLYPH_HEADER

void main()
{
  vec4 fColor = texture2D(textureMap, textureLookup + gl_PointCoord * textureScale);
#include ANAGLYPH_BODY
}

