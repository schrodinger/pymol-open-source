
// cylinder imposter fragment shader

uniform bool lighting_enabled;

uniform float fog_enabled;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

uniform bool bg_gradient;
uniform float inv_height;
uniform float ortho;
uniform float no_flat_caps;
uniform bool filter_front_facing;
uniform bool two_sided_lighting_enabled;
uniform int light_count;
uniform float shininess;
uniform float shininess_0;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;
uniform float half_bond;

#include ANAGLYPH_HEADER

//varying vec3 point; // surface point
//varying vec3 axis; // cylinder axis
//varying vec3 base; // cylinder base
//varying vec3 end; // cylinder end
//varying vec3 U; // cylinder base plane coordinates
//varying vec3 V;
//varying float radius; // radius
//varying float cap; // should we draw the endcap
//varying float inv_sqr_height;

varying vec4 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_3 ;
varying vec4 packed_data_4 ;
varying vec4 packed_data_5 ;

// point -> surface_point b/c preprocessor replaces _point
#define surface_point ( packed_data_0.xyz )
#define axis ( packed_data_1.xyz )
#define base ( packed_data_2.xyz )
// end -> end_cyl
#define end_cyl packed_data_3.xyz
#define U ( packed_data_4.xyz )
#define V ( packed_data_5.xyz )
#define radius ( packed_data_3.w )
#define cap ( packed_data_4.w )
#define inv_sqr_height ( packed_data_5.w )

varying vec4 color1;
varying vec4 color2;

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

#include ComputeColorForLight

void main(void)
{
#ifndef cylinder_shader_ff_workaround
    // cull back face - otherwise we are drawing all pixels twice
    // this change gives roughly 2x speedup
    if (filter_front_facing && !gl_FrontFacing) 
      discard; 
#endif

    vec3 ray_target = surface_point;
    vec3 ray_origin = vec3(0.0);
    vec3 ray_direction = mix(normalize(ray_origin - ray_target), vec3(0.0, 0.0, 1.0), ortho);

    // basis is local system of coordinates for the cylinder
    mat3 basis = mat3(U, V, axis);

    vec3 diff = ray_target - 0.5 * (base + end_cyl);
    vec3 P = diff * basis;

    // angle (cos) between cylinder cylinder_axis and ray direction
    float dz = dot(axis, ray_direction);

    float radius2 = radius*radius;

    // calculate distance to the cylinder from ray origin
    vec3 D = vec3(dot(U, ray_direction),
                  dot(V, ray_direction),
                  dz);
    float a0 = P.x*P.x + P.y*P.y - radius2;
    float a1 = P.x*D.x + P.y*D.y;
    float a2 = D.x*D.x + D.y*D.y;
    // calculate a dicriminant of the above quadratic equation
    float d = a1*a1 - a0*a2;
    if (d < 0.0)
        // outside of the cylinder
        discard;

    float dist = (-a1 + sqrt(d))/a2;

    // point of intersection on cylinder surface
    vec3 new_point = ray_target + dist * ray_direction;

    vec3 tmp_point = new_point - base;
    vec3 normal = normalize(tmp_point - axis * dot(tmp_point, axis));

    float ratio = dot(new_point-base, vec3(end_cyl-base)) * inv_sqr_height;

    ray_origin = mix(ray_origin, surface_point, ortho);

    /* cap :  4 bits : 1st - frontcap
                       2nd - endcap
                       3rd - frontcapround
                       4th - endcapround
		       5th - smooth_half_bonds
     */
    float fcap = cap;
    float frontcap = 0.0, frontcapround = 0.0;
    float endcap = 0.0, endcapround = 0.0;
    float smooth_half_bonds = 0.0;

    frontcap = mod(fcap, 2.0);
    fcap = ( fcap - frontcap ) / 2.0;
    endcap = mod(fcap, 2.0);
    fcap = ( fcap - endcap ) / 2.0;
    frontcapround = floor((mod(fcap, 2.0) + no_flat_caps) / 2.0);
    fcap = ( fcap - frontcapround ) / 2.0;
    endcapround = floor((mod(fcap, 2.0) + no_flat_caps) / 2.0);
    fcap = ( fcap - endcapround ) / 2.0;
    
    smooth_half_bonds = floor((mod(fcap, 2.0) + no_flat_caps) / 2.0);

    vec4 color;

    float dp = clamp(-half_bond*new_point.z*inv_height, 0., .5) * smooth_half_bonds;
    color = mix(color1, color2, smoothstep(.5 - dp, .5 + dp, ratio));

    // test front cap
    float cap_test = dot((new_point - base), axis);

    // to calculate caps, simply check the angle between
    // the point of intersection - cylinder end vector
    // and a cap plane normal (which is the cylinder cylinder_axis)
    // if the angle < 0, the point is outside of cylinder
    // test front cap

    // flat
    if (frontcapround < 0.5 && cap_test < 0.0) {
      // ray-plane intersection
      color = color1;
      float dNV = dot(-axis, ray_direction);
      if (dNV < 0.0) discard;
      float near = dot(-axis, (base)) / dNV;
      new_point = ray_direction * near + ray_origin;
      // within the cap radius?
      if (dot(new_point - base, new_point-base) > radius2) discard;
      normal = -axis;
    }

    // round
    if (frontcapround > 0.5 && cap_test < 0.0) {
      if ( frontcap < 0.5)
        discard;
      color = color1;
      vec3 sphere_direction = mix(base, ray_origin - base, ortho);
      float b = dot(sphere_direction, ray_direction);
      float pos = b*b + radius2 -dot(sphere_direction, sphere_direction);
      if ( pos < 0.0)
        discard;
      float near = mix(b + sqrt(pos), sqrt(pos) - b, ortho);
      new_point = near * ray_direction + ray_origin;     
      normal = normalize( new_point - base ); 
    }
    // test end cap

    cap_test = dot((new_point - end_cyl), axis);

    // flat
    if (endcapround < 0.5 && cap_test > 0.0) {
      // ray-plane intersection
      color = color2;
      float dNV = dot(axis, ray_direction);
      if (dNV < 0.0) discard;
      float near = dot(axis, end_cyl) / dNV;
      new_point = ray_direction * near + ray_origin;
      // within the cap radius?
      if (dot(new_point - end_cyl, new_point-base) > radius2) discard;
      normal = axis;
    }

    // round

    if (endcapround > 0.5 && cap_test > 0.0) {
      if ( endcap < 0.5)
        discard;
      color = color2;
      vec3 sphere_direction = mix(end_cyl, ray_origin - end_cyl, ortho);
      float b = dot(sphere_direction, ray_direction);
      float pos = b*b + radius2 -dot(sphere_direction, sphere_direction);
      if ( pos < 0.0)
        discard;
      float near = mix(b + sqrt(pos), sqrt(pos) - b, ortho);
      new_point = near * ray_direction + ray_origin;
      normal = normalize( new_point - end_cyl );
    }

    vec2 clipZW = new_point.z * gl_ProjectionMatrix[2].zw +
        gl_ProjectionMatrix[3].zw;

    float depth = 0.5 + 0.5 * clipZW.x / clipZW.y;

    // this is a workaround necessary for Mac
    // otherwise the modified fragment won't clip properly

    if (depth <= 0.0)
      discard;

    if (depth >= 1.0)
      discard;

    gl_FragDepth = depth;

  vec4 final_color;
  int i;
  float NdotL, NdotH;

  final_color = (gl_LightModel.ambient) * color;

#include CallComputeColorForLight

  float cfog = clamp((g_Fog_end + new_point.z) * g_Fog_scale, 0.0, 1.0);

  cfog = mix(1.0, clamp(cfog, 0.0, 1.0), fog_enabled);

  vec4 fogColor = ComputeFogColor();

  final_color.rgb = mix(fogColor.xyz, final_color.rgb, cfog);

  vec4 fColor = vec4(final_color.rgb, color.a);

#include ANAGLYPH_BODY
}

