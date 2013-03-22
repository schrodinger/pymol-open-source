/* ios recent */
uniform sampler2D textureMap;
varying vec2 textureLookup ;
varying vec4 pickcolor ;
uniform float isPicking;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;
uniform float fog_enabled;
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
  float cfog = clamp(mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled), 0.0, 1.0); // not sure why extra clamp is needed on ipad, not simulator
  vec4 finalColor = texture2D(textureMap, textureLookup);
  vec4 fogColor = ComputeFogColor();
  finalColor.xyz = mix(fogColor.xyz, finalColor.xyz, finalColor.a);
  vec4 fColor = vec4(mix(fogColor.xyz, finalColor.xyz, cfog), finalColor.a);
#ifdef PYMOL_IOS
  vec4 npColor = (1. - isPicking) * fColor;
  if (npColor.a < .1)
     discard;
  gl_FragColor = npColor + isPicking * pickcolor;
#else
  gl_FragColor = (1. - isPicking) * fColor + isPicking * pickcolor;
#endif
}

