#ifdef ANAGLYPH
  fColor.rgb = matL * pow(fColor.rgb, vec3(gamma));
#endif
gl_FragColor = fColor;
