/*
 * Lighting fragment header
 */

uniform float shininess;
uniform float shininess_0;
uniform float spec_value;
uniform float spec_value_0;

#ifdef precomputed_lighting
uniform samplerCube lightingTex;
#else
uniform struct {
  vec4 ambient;
  vec4 diffuse;
  vec4 specular;
  vec4 position;
} g_LightSource[8];

uniform struct {
  vec4 ambient;
} g_LightModel;
#endif

#ifdef ray_transparency_oblique
uniform float trans_oblique;
uniform float oblique_power;
#endif

/*
 * Call this function for every light
 *
 * Return: (color intensity, specular (white) intensity)
 */
vec2 ComputeLighting(vec3 normal,   // surface normal
                     vec3 L,        // gl_LightSource[i].position.xyz
                     float diffuse, // gl_LightSource[i].diffuse.r
                     float spec,    // specular intensity
                     float shine)   // specular exponent
{
  // light direction (normalized)
  L = normalize(L);

  // cosine of angle between normal and light
  float NdotL = dot(normal, L);

  // normals that don't point away from the light
  if (NdotL > 0.0) {

  // diffuse
  diffuse *= NdotL;

  // specular
  vec3 H = normalize(L + vec3(0., 0., 1.));
  float NdotH = max(dot(normal, H), 0.0);
  spec *= pow(NdotH, shine);

  return vec2(diffuse, spec);

  } else {
    return vec2(0.0);
  }
}

/*
 * Apply lighting from all lights
 *
 * Return: lighted color
 */
vec4 ApplyLighting(vec4 color, vec3 normal) {
#ifdef precomputed_lighting
  vec2 lighting = textureCube(lightingTex, normal).ra;
#else
  vec2 lighting = vec2(g_LightModel.ambient.r, 0.0);

  // add to lighting for every light
#include CallComputeColorForLight
#endif

  color.rgb *= min(lighting.x, 1.0);
  color.rgb += lighting.y;

#ifdef ray_transparency_oblique
  // see Ray.cpp
  float oblique_factor = abs(normal.z);
  oblique_factor = mix(
      (0.5 + 0.5 * (1.0 -
                    pow((1.0 - oblique_factor) * 2.0, oblique_power))),
      (      0.5 * (
                    pow((      oblique_factor) * 2.0, oblique_power))),
      step(oblique_factor, 0.5));
  float trans = (1. - color.a) * (trans_oblique * oblique_factor + (1. - trans_oblique));
  color.a = 1. - clamp(trans, 0.0, 1.0);
#endif

  return color;
}
