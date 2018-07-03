/*
 * Background and fog fragment header
 */

uniform sampler2D bgTextureMap;

uniform vec3 bgSolidColor;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;

vec3 ComputeBgColor() {
#ifdef bg_image_mode_solid
  return bgSolidColor.rgb;
#else

  float covered = 1.0;

#ifdef bg_image_mode_1_or_3
  // mode 1 (or 3): centered (and repeated)
  vec2 bgLookup = 0.5 + (bgTextureLookup - 0.5) / viewImageSize;
#ifndef bg_image_mode_2_or_3
  // mode 1: centered
  vec2 within01 = step(0., bgLookup) * step(-1., -bgLookup);
  covered = within01.x * within01.y;
#endif
#else
#ifdef bg_image_mode_2_or_3
  // mode 2: tiled
  vec2 bgLookup = mod(bgTextureLookup, tiledSize) / tiledSize;
  bgLookup = tileSize * floor(bgLookup / tileSize);
#else
  // mode 0: stretched
  vec2 bgLookup = bgTextureLookup;
#endif
#endif

  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  return mix(bgSolidColor.rgb, bgColor.rgb, min(covered, bgColor.a));
#endif
}

vec4 ApplyFog(vec4 color, float fog) {
#ifdef depth_cue
  if (isPicking || fog >= 1.)
    return color;

  return vec4(mix(ComputeBgColor(), color.rgb, fog), color.a);
#else
  return color;
#endif
}
