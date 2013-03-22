/* ios recent */
uniform sampler2D textureMap;

varying vec4 backgroundColor;
varying vec2 textureLookup ;

#include ANAGLYPH_HEADER

void main()
{
  vec4 fColor = texture2D(textureMap, textureLookup);
#ifdef PYMOL_IOS
  if (isLabel > .5 && fColor.a < .1)
     discard;
#endif
  vec4 bColor = backgroundColor * (1.-fColor.a);
  fColor = vec4(bColor.rgb + fColor.rgb * fColor.a, bColor.a + fColor.a );

#include ANAGLYPH_BODY
}

