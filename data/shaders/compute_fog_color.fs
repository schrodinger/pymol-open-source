vec4 ComputeFogColor(){
#ifdef bg_image_mode_stretched
  vec2 bgLookup = pixelSize * floor(bgTextureLookup / pixelSize);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  return vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, 1.);
#endif
#ifdef bg_image_mode_solid
  return vec4(fogSolidColor.rgb, 1.);	
#endif
}
