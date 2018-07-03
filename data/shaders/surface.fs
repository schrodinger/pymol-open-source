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
uniform vec4 interior_color;
uniform bool use_interior_color;

#include anaglyph_header.fs
#include compute_fog_color.fs
#include compute_color_for_light.fs

void main()
{
  if (isPicking){
    gl_FragColor = COLOR; // no lighting
  } else {
    // WARNING: !gl_FrontFacing does not work on Intel (bug)
    bool is_interior = gl_FrontFacing ? false : true;
    vec3 normal = NORMAL;
    if (use_interior_color && is_interior) {
      // ray_interior_color
      gl_FragColor = vec4(interior_color.rgb, COLOR.a);
    } else {
      // back faces
      if (is_interior) {
        if (two_sided_lighting_enabled) {
          normal = -normal;
        } else {
          // disable all except ambient light (intended? ray tracing differs!)
          normal = vec3(0.);
        }
      }

      vec4 color = ApplyColorEffects(COLOR, gl_FragCoord.z);

      // lights
      color = ApplyLighting(color, normal);

      // fog
      gl_FragColor = ApplyFog(color, fog);
    }
    PostLightingEffects();
  }
}

