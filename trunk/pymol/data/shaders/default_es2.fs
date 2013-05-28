
varying vec3 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_3 ;
varying vec4 packed_data_4 ;

//varying vec3 N;
#define COLOR packed_data_3
#define NORMAL normalize(packed_data_0.xyz)

uniform float fog_enabled;
varying float fog;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

varying vec2 bgTextureLookup;

uniform bool lighting_enabled;
uniform bool two_sided_lighting_enabled;
uniform int light_count;
uniform vec4 interior_color;
uniform float interior_color_threshold;
uniform float shininess;
uniform float shininess_0;
uniform bool use_interior_color_threshold;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;

#include ANAGLYPH_HEADER

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

#include ComputeColorForLight

void main()
{
  vec4 final_color = vec4(0.);

  if (lighting_enabled){
    bool is_interior = false;
    if (use_interior_color_threshold || two_sided_lighting_enabled){
      vec3 viewV = vec3(0.,0.,-1.);
      float dotp = dot(NORMAL, viewV);
      is_interior = ( dotp > interior_color_threshold );
    }
    if (!two_sided_lighting_enabled && is_interior){
      final_color = interior_color;
    } else {
      final_color = (gl_LightModel.ambient) * COLOR;

#include CallComputeColorForLight
    }
  } else {
    final_color = COLOR;
  }

  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  vec4 fogColor = ComputeFogColor();
  vec4 fColor = vec4(mix(vec3(fogColor), final_color.rgb, cfog), COLOR.a);

#include ANAGLYPH_BODY
}

