/* ios recent */
uniform sampler2D textureMap;
varying vec2 textureLookup ;

void main()
{
  vec4 fColor = texture2D(textureMap, textureLookup);
#ifdef PYMOL_IOS
  if (fColor.a < .1)
     discard;
#endif
  gl_FragColor = fColor;
}

