
uniform bool isPicking;

#ifdef ANAGLYPH
uniform mat3 matL;
uniform float gamma;
#endif

#ifdef NO_ORDER_TRANSP

#ifdef ONE_DRAW_BUFFER
uniform float which_pass;  // 0. - first, 1. - second
#endif

#define gl_FragColor gl_FragData[0]

float get_weight(float depth, float alpha){
      return alpha * max(1e-2, 3e3 * pow(1. - depth, 3.));  // eq #10
}

/*
 * Set data for order independant transparency
 *
 * With multiple buffers:
 * gl_FragData[0]: weighted color
 * gl_FragData[1]: weight
 *
 * With one buffer and multi-pass:
 * gl_FragData[0]: weighted color in pass 1
 * gl_FragData[0]: weight         in pass 2
 *
 */
void setFragDataForNoOrderTransp(float depth) {
  float wi = get_weight(depth, gl_FragColor.a);

  gl_FragData[0] = vec4(gl_FragColor.rgb * wi, gl_FragColor.a);

#ifdef ONE_DRAW_BUFFER
  // for using one draw buffer, switch on which_pass between color and weights
  gl_FragData[0] = mix(gl_FragData[0], vec4(wi, wi, 0., 1.), which_pass);
#else
#ifdef ANAGLYPH
  // need to set two components because anaglyph masks out colors for each eye
  gl_FragData[1].rg = vec2(wi);
#else
  gl_FragData[1].r = wi;
#endif // ANAGLYPH
#endif // ONE_DRAW_BUFFER
}

#endif // NO_ORDER_TRANSP

/*
 * Apply post-lighting effects directly to output variables like gl_FragColor
 */
void PostLightingEffectsNoOIT() {
  if (isPicking)
    return;
#ifdef ray_trace_mode_3
  // cel shading
  gl_FragColor.rgb = floor(gl_FragColor.rgb * 3.999) / 3.0;
#endif

#ifdef ANAGLYPH
  gl_FragColor.rgb = matL * pow(gl_FragColor.rgb, vec3(gamma));
#endif
}

vec3 ComputeChromadepth(float depth);

/*
 * ... including OIT
 */
void PostLightingEffects(float depth) {
  if (isPicking)
    return;
  PostLightingEffectsNoOIT();

#ifdef chromadepth_postlighting
  vec3 luma_coeff = vec3(.30, .59, .11);
  gl_FragColor.rgb = ComputeChromadepth(depth) * dot(luma_coeff, gl_FragColor.rgb);
#endif

#ifdef transparency_mode_3
#ifndef NO_ORDER_TRANSP
  // if transparent, discard in opaque pass
  if (gl_FragColor.a < .95)
    discard;
#else
  // if almost opaque, discard in transparent pass
  if (gl_FragColor.a > .95)
    discard;

  setFragDataForNoOrderTransp(depth);
#endif // NO_ORDER_TRANSP
#endif // transparency_mode_3
}

/*
 * Default argument overload
 */
void PostLightingEffects() {
  PostLightingEffects(gl_FragCoord.z);
}

#ifdef chromadepth
#ifndef ortho
/*
 * "unproject" the depth value if not in orthoscopic view, to get
 * a depth gradient which is linear in eye space.
 */
uniform vec2 clippingplanes;
float orthoDepth(float depth) {
  float front = clippingplanes.x;
  float back = clippingplanes.y;
  float front_inv = 1. / front;
  float back_inv  = 1. / back;

  // from proj normalized dev coords to eye space
  float depth_eye = 1. / (depth * (back_inv - front_inv) + front_inv);

  // from eye space to ortho normalized dev coords
  return (depth_eye - front) / (back - front);
}
#endif

/*
 * Chromadepth - Calculates a color by depth.
 * Uses a rainbow from red (near) to blue (far).
 *
 * http://chromatek.com/art-and-design/design-guide/
 */
vec3 ComputeChromadepth(float depth) {
#ifndef ortho
  depth = orthoDepth(depth);
#endif

  const float margin = 1. / 4.;
  float hue6 = 6. * 240. / 360. *
    clamp(depth * (1. + 2. * margin) - margin, 0., 1.);

  return clamp(vec3(
        abs(hue6 - 3.) - 1.,
        2. - abs(hue6 - 2.),
        2. - abs(hue6 - 4.)),
      0., 1.);
}
#endif

/*
 * Applies pre-lighting and pre-fog color alteration effects.
 */
vec4 ApplyColorEffects(vec4 color, float depth) {
#ifdef chromadepth
#ifndef chromadepth_postlighting
  if (!isPicking)
    color.rgb = ComputeChromadepth(depth);
#endif
#endif
  return color;
}
