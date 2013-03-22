// Sphere impostor fragment shader

uniform bool lighting_enabled;

uniform float ortho;

uniform float fog_enabled;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

uniform bool bg_gradient;
uniform float inv_height;
uniform int light_count;
uniform float shininess;
uniform float shininess_0;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;

#include ANAGLYPH_HEADER

varying vec4 COLOR;
varying vec3 sphere_center;
varying float radius2;
varying vec3 point;

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
    vec3 ray_origin = mix(vec3(0.,0.,0.), point, ortho);

    vec3 ray_direction = mix(normalize(point), vec3(0., 0., 1.), ortho);

    vec3 sphere_direction = mix(sphere_center, ray_origin - sphere_center, ortho);

    // Calculate sphere-ray intersection
    float b = dot(sphere_direction, ray_direction);

    float position = b * b + radius2 - dot(sphere_direction, sphere_direction);

    // Check if the ray missed the sphere
    if (position < 0.0)
       discard;

    // Calculate nearest point of intersection
    float nearest = mix(b - sqrt(position), sqrt(position) - b, ortho);

    // Calculate intersection point on the sphere surface.  The ray
    // origin is at the quad (center point), so we need to project
    // back towards the user to get the front face.
    vec3 ipoint = nearest * ray_direction + ray_origin;

    // Calculate normal at the intersection point
    vec3 N = normalize(ipoint - sphere_center);

    // Calculate depth in clipping space 
    vec2 clipZW = ipoint.z * gl_ProjectionMatrix[2].zw +
        gl_ProjectionMatrix[3].zw;

    float depth = 0.5 + 0.5 * clipZW.x / clipZW.y;

    // this is a workaround necessary for Mac
    // otherwise the modified fragment won't clip properly

/*
    float isDiscarded = step(.5, step(depth, 0.) + step(1.-depth, 0.));
    if (isDiscarded > 0.0)
      discard;
*/
    if (depth <= 0.0)
      discard;

    if (depth >= 1.0)
      discard;

    gl_FragDepth = depth;

    vec4 color;

    vec3 L0 = normalize(vec3(gl_LightSource[0].position) - ipoint);

    float NdotL = max(dot(N, L0), 0.0);
    float NdotH;

    vec4 final_color = (gl_LightModel.ambient) * COLOR;

#include CallComputeColorForLight

/*
    int i;
    for (i=0; i<light_count;i++){
      vec3 L = normalize(gl_LightSource[i].position.xyz);
      vec3 H = normalize(gl_LightSource[i].halfVector.xyz);
      float spec = 0., shine = 0.;
      if (i==0){
        spec = spec_value_0;
        shine = shininess_0;
      } else if (spec_count >= i){
        spec = spec_value;
        shine = shininess;
      }
      final_color += gl_LightSource[i].ambient * COLOR;
      NdotL = dot(N, L);
      if (NdotL > 0.0) {
        final_color += gl_LightSource[i].diffuse * NdotL * COLOR;
        NdotH = max(dot(N, H), 0.0);
        final_color += spec * pow(NdotH, shine);
      }
    }
*/
    float fogv = (g_Fog_end + ipoint.z) * g_Fog_scale;
    float cfog = clamp(fogv, 0.0, 1.0);

    cfog = mix(1.0, clamp(cfog, 0.0, 1.0), fog_enabled);

    vec4 fogColor = ComputeFogColor();

    final_color.rgb = mix(fogColor.rgb, final_color.rgb, cfog);

    vec4 fColor = vec4(final_color.rgb, COLOR.a);

#include ANAGLYPH_BODY
}
