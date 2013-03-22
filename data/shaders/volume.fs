uniform sampler3D volumeTex;
uniform sampler1D colorTex;

uniform float fog_enabled;
uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;
varying float fog;

uniform float g_Fog_end;
uniform float g_Fog_scale;
varying vec2 bgTextureLookup;

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
  vec4 color = texture1D(colorTex, texture3D(volumeTex, gl_TexCoord[0].xyz).x);
  if (color.a == 0.0) discard;
  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  vec4 fogColor = ComputeFogColor();
  gl_FragColor = vec4(vec3(mix(fogColor.rgb, color.rgb, cfog)), color.a);
}

