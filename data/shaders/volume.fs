uniform sampler3D volumeTex;
uniform sampler2D colorTex2D;
uniform sampler1D colorTex1D;
uniform float volumeScale;
uniform float volumeBias;
uniform float sliceDist;
uniform mat4 TexMatrix;
uniform vec3 eyeposM;
uniform vec3 vDirM;
uniform sampler3D carvemask;
uniform bool carvemaskFlag;

varying vec3 vertexM;
varying float fog;
varying vec2 bgTextureLookup;

#include anaglyph_header.fs
#include compute_fog_color.fs

bool iscarvemasked(vec3 t) {
  return carvemaskFlag && texture3D(carvemask, t).r > 0.5;
}

void main()
{
#ifdef volume_mode
#else // volume_mode

  if (iscarvemasked(gl_TexCoord[0].xyz))
    discard;

  float v = texture3D(volumeTex, gl_TexCoord[0].xyz).r;
  v = v * volumeScale + volumeBias;
  if (v < 0. || v > 1.) discard;
  vec4 color = texture1D(colorTex1D, v);
#endif // volume_mode

  if (color.a == 0.0)
    discard;

  color = ApplyColorEffects(color, gl_FragCoord.z);

  gl_FragColor = ApplyFog(color, fog);

  PostLightingEffects();
}

