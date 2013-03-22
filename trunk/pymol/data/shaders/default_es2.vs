attribute vec4 a_Vertex;
attribute vec3 a_Normal;
attribute vec4 a_Color;
attribute float a_Accessibility; /* This is for ambient occlusion, 1.0 by default */

varying vec3 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_3 ;
varying vec4 packed_data_4 ;

//varying vec3 N;
#define NORMAL  packed_data_0
#define COLOR packed_data_3
//#define fog packed_data_1.w
#define E 2.718281828459045
/* this doesn't work for setting fog_color, need to set the values independently */
//#define fog_color ( packed_data_1.x, packed_data_1.y, packed_data_1.z )

varying float fog;
varying vec2 bgTextureLookup;
uniform vec2 t2PixelSize;

uniform bool lighting_enabled;

uniform bool bg_gradient;
uniform float ambient_occlusion_scale;
uniform int accessibility_mode;

void main()
{
  NORMAL = normalize(gl_NormalMatrix * a_Normal);
  vec3 eye_pos = vec3(gl_ModelViewMatrix * a_Vertex);
  if (accessibility_mode == 1){
    COLOR = vec4(clamp(a_Color.xyz * (1.-(ambient_occlusion_scale*a_Accessibility)), 0., 1.), a_Color.w);
  } else if (accessibility_mode == 2){
    COLOR = vec4(a_Color.xyz * cos(90.*radians(clamp(ambient_occlusion_scale*a_Accessibility, 0., 1.))), a_Color.w);
  } else {
    COLOR = vec4(clamp(a_Color.xyz * (1. / (1. + pow(E, 0.5 * ( (ambient_occlusion_scale * a_Accessibility) - 10.)))), 0., 1.), a_Color.w);
  }
  
  // This was breaking fog on ATI/Linux
  //  gl_FogFragCoord = -eye_pos.z;
  //  fog = (gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale;
  fog = (gl_Fog.end - abs(eye_pos.z)) * gl_Fog.scale;
  gl_Position = gl_ModelViewProjectionMatrix * a_Vertex;
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
