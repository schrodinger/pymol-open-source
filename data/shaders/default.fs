#include webgl_header.fs

varying vec3 packed_data_0 ;
varying float fog;
varying vec2 bgTextureLookup;

#ifdef PURE_OPENGL_ES_2
varying vec4 COLOR;
#else
// using the built-in allows to use glShadeModel
#define COLOR gl_Color
#endif

#define NORMAL normalize(packed_data_0.xyz)

uniform bool lighting_enabled;
uniform bool two_sided_lighting_enabled;

#include anaglyph_header.fs
#include compute_fog_color.fs
#include compute_color_for_light.fs

void main()
{
  vec4 color = ApplyColorEffects(COLOR, gl_FragCoord.z);

  if (!isPicking && lighting_enabled){
    vec3 normal = NORMAL;
    if (two_sided_lighting_enabled){
       bool ff = gl_FrontFacing ? true : false;
       if (!ff)
          normal = -normal;
    }
    color = ApplyLighting(color, normal);
  }

  gl_FragColor = ApplyFog(color, fog);

  PostLightingEffects();
}

