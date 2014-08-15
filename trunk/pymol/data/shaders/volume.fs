uniform sampler3D volumeTex;
uniform sampler1D colorTex;
uniform float volumeScale;
uniform float volumeBias;
uniform sampler3D carvemask;
uniform bool carvemaskFlag;

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

bool iscarvemasked(vec3 t) {
  return carvemaskFlag && texture3D(carvemask, t).r > 0.5;
}

void main()
{
  if (iscarvemasked(gl_TexCoord[0].xyz))
    discard;

  float v = texture3D(volumeTex, gl_TexCoord[0].xyz).r;
  v = v * volumeScale + volumeBias;
  if (v < 0. || v > 1.) discard;
  vec4 color = texture1D(colorTex, v);

  if (color.a == 0.0) discard;
  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  vec4 fogColor = ComputeFogColor();
  gl_FragColor = vec4(vec3(mix(fogColor.rgb, color.rgb, cfog)), color.a);
}

