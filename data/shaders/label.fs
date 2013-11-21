/* ios recent */
uniform sampler2D textureMap;
varying vec2 textureLookup ;
varying vec4 pickcolor ;
uniform float isPicking;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;
varying float fog;
varying vec3 normalizedViewCoordinate;
#define bgTextureLookup normalizedViewCoordinate.xy

uniform float isStretched;
uniform float isCentered;
uniform float isCenteredOrRepeated;
uniform float isTiled;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;

#include ComputeFogColor

void main()
{
  vec4 fColor = texture2D(textureMap, textureLookup);
  if (fog < 1.0) {
      vec4 fogColor = ComputeFogColor();
      fColor.rgb = mix(fogColor.rgb, fColor.rgb, fog);
  }
#ifdef PYMOL_IOS
  vec4 npColor = (1. - isPicking) * fColor;
  if (npColor.a < .1)
     discard;
  gl_FragColor = npColor + isPicking * pickcolor;
#else
  gl_FragColor = (1. - isPicking) * fColor + isPicking * pickcolor;
#endif
}

