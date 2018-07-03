#include webgl_header.vs

uniform bool isPicking;

attribute vec4 a_Vertex;
attribute vec3 a_Normal;
attribute vec4 a_Color;
attribute float a_Accessibility; /* This is for ambient occlusion, 1.0 by default */

varying vec3 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_4 ;

//varying vec3 N;
#define NORMAL  packed_data_0

#ifdef PURE_OPENGL_ES_2
varying vec4 COLOR;
#else
// using the built-in allows to use glShadeModel
#define COLOR gl_FrontColor
#endif

//#define fog packed_data_1.w
#define E 2.718281828459045
/* this doesn't work for setting fog_color, need to set the values independently */
//#define fog_color ( packed_data_1.x, packed_data_1.y, packed_data_1.z )

varying float fog;
varying vec2 bgTextureLookup;

uniform bool lighting_enabled;

uniform float ambient_occlusion_scale;
uniform int accessibility_mode;
uniform float accessibility_mode_on;

void main()
{
  NORMAL = normalize(g_NormalMatrix * a_Normal);
  vec3 eye_pos = vec3(g_ModelViewMatrix * a_Vertex);

#ifdef PYMOL_WEBGL_IOS
  gl_PointSize = g_PointSize;
#endif

  if (isPicking){
    COLOR = a_Color; // no lighting
  } else {
    vec4 COLORa;
    if (accessibility_mode == 1){
      COLORa = vec4(clamp(a_Color.xyz * (1.-(ambient_occlusion_scale*a_Accessibility)), 0., 1.), a_Color.w);
    } else if (accessibility_mode == 2){
      COLORa = vec4(a_Color.xyz * cos(90.*radians(clamp(ambient_occlusion_scale*a_Accessibility, 0., 1.))), a_Color.w);
    } else {
      COLORa = vec4(clamp(a_Color.xyz * (1. / (1. + pow(E, 0.5 * ( (ambient_occlusion_scale * a_Accessibility) - 10.)))), 0., 1.), a_Color.w);
    }
    COLOR = mix(a_Color, COLORa, accessibility_mode_on);
    fog = (g_Fog_end + eye_pos.z) * g_Fog_scale;
  }
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
