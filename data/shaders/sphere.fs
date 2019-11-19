#ifdef PYMOL_WEBGL_IOS
#extension GL_EXT_frag_depth : require
#endif

#include webgl_header.fs

// Sphere impostor fragment shader

varying vec4 COLOR;
varying vec3 sphere_center;
varying float radius2;
varying vec3 point;

varying vec2 bgTextureLookup;

uniform bool lighting_enabled;
uniform mat4 g_ProjectionMatrix;
uniform float g_Fog_end;
uniform float g_Fog_scale;

#include anaglyph_header.fs
#include compute_fog_color.fs
#include compute_color_for_light.fs

void main(void)
{
#ifdef ortho
    vec3 ray_origin = point;
    vec3 ray_direction = vec3(0., 0., -1.);
    vec3 sphere_direction = ray_origin - sphere_center;
#else
    vec3 ray_origin = vec3(0., 0., 0.);
    vec3 ray_direction = normalize(point);
    vec3 sphere_direction = sphere_center;
#endif

    // Calculate sphere-ray intersection
    float b = dot(sphere_direction, ray_direction);

    float position = b * b + radius2 - dot(sphere_direction, sphere_direction);

    // Check if the ray missed the sphere
    if (position < 0.0)
       discard;

    // Calculate nearest point of intersection
    float nearest = b - sqrt(position);

    // Calculate intersection point on the sphere surface.  The ray
    // origin is at the quad (center point), so we need to project
    // back towards the user to get the front face.
    vec3 ipoint = nearest * ray_direction + ray_origin;

    // Calculate normal at the intersection point
    vec3 normal = normalize(ipoint - sphere_center);

    // Calculate depth in clipping space 
    vec2 clipZW = ipoint.z * g_ProjectionMatrix[2].zw +
        g_ProjectionMatrix[3].zw;

    float depth = 0.5 + 0.5 * clipZW.x / clipZW.y;

    // this is a workaround necessary for Mac
    // otherwise the modified fragment won't clip properly

/*
    float isDiscarded = step(.5, step(depth, 0.) + step(1.-depth, 0.));
    if (isDiscarded > 0.0)
      discard;
*/
    if (depth <= 0.0 || depth >= 1.0)
      discard;

    gl_FragDepth = depth;

    if (!isPicking){
      if (lighting_enabled){
        vec4 color = ApplyColorEffects(COLOR, depth);
        color = ApplyLighting(color, normal);
        float fogv = (g_Fog_end + ipoint.z) * g_Fog_scale;
        gl_FragColor = ApplyFog(color, fogv);
      } else {
        gl_FragColor = COLOR;
      }
      PostLightingEffects(depth);
    } else {
      if (COLOR.a == 0.0) {
        // cPickableThrough
        discard;
      }
      gl_FragColor = COLOR;
    }

}
