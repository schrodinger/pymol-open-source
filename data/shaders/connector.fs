#include webgl_header.fs

varying vec3 NORMAL ;
#ifdef use_geometry_shaders
varying vec4 COLOR ;
#else
varying vec4 COLORIn ;
#define COLOR COLORIn
#endif

varying float fog;
varying vec2 bgTextureLookup;

#ifdef use_geometry_shaders
varying float lineEdge;
varying float aaCutoff;
#endif

#include anaglyph_header.fs
#include compute_fog_color.fs

float mysmoothstep(float edge0, float edge1, float x){
  float rets = step(edge0, edge1) * step(edge1, edge0);
  return rets * step(edge0, x) + (1.-rets) * smoothstep(edge0, edge1, x);
}

void main()
{
  float alpha = COLOR.a;

#ifdef use_geometry_shaders
#ifdef line_smooth
  alpha *= mysmoothstep(0., aaCutoff, 1. - abs(lineEdge));
#endif
#endif

  vec4 color = ApplyColorEffects(COLOR, gl_FragCoord.z);

  gl_FragColor = ApplyFog(vec4(color.rgb, alpha), fog);

  PostLightingEffects();
}

