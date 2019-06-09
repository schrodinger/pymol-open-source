#include webgl_header.fs

#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

#ifdef line_smooth
uniform float line_width;
varying float centerdist;
#endif

varying vec4 COLOR;
varying vec4 COLOR2;
varying float a_INTERPOLATE;
varying float fog;
varying vec2 bgTextureLookup;
varying float whichEnd;

#include anaglyph_header.fs
#include compute_fog_color.fs

void main(void) {
  float which = mix(step(0.5, whichEnd), whichEnd, a_INTERPOLATE);
  vec4 BCOLOR = mix(COLOR, COLOR2, which);
  vec4 color = ApplyColorEffects(BCOLOR, gl_FragCoord.z);

#ifdef line_smooth
  const float margin = 1.5;
  if (!isPicking)
    color.a *= (1. - max(0.,
                (abs(centerdist) - line_width + margin) / margin));
#endif

  gl_FragColor = ApplyFog(color, fog);

  PostLightingEffectsNoOIT();
}
