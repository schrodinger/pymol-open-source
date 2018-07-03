#include webgl_header.fs

//#ifndef use_geometry_shaders

#define bgTextureLookup a_bgTextureLookup

varying vec3 a_NORMAL ;
varying vec4 a_COLOR ;
varying vec4 a_COLOR2 ;
varying vec4 a_COLOR_INTERP ;
varying float a_INTERPOLATE;
varying float a_LINE_POSITION;
varying float a_fog;
varying vec2 a_bgTextureLookup;
 
uniform bool lighting_enabled;
uniform bool two_sided_lighting_enabled;

#include anaglyph_header.fs
#include compute_fog_color.fs
#include compute_color_for_light.fs

void main()
{
  float whichColor = step(.5, a_LINE_POSITION);
  vec4 color_step = whichColor * a_COLOR2 / a_LINE_POSITION + (1.-whichColor) * a_COLOR / (1.-a_LINE_POSITION);
  vec4 icolor = a_INTERPOLATE * a_COLOR_INTERP + (1.-a_INTERPOLATE) * color_step;
  vec4 color = ApplyColorEffects(icolor, gl_FragCoord.z);
  if (!isPicking && lighting_enabled){
    vec3 normal = a_NORMAL;
    if (two_sided_lighting_enabled){
       bool ff = gl_FrontFacing ? true : false;
       if (!ff)
          normal = -normal;
    }
    color = ApplyLighting(color, normal);
  }

  gl_FragColor = ApplyFog(color, a_fog);

  PostLightingEffects();
}

