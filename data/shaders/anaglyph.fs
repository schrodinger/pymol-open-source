#ifdef ANAGLYPH
  gl_FragColor = stereo_flag_l * vec4(matL * pow(fColor.rgb,vec3(gamma,gamma,gamma)), fColor.a)
               + stereo_flag_r * vec4(matR * pow(fColor.rgb,vec3(gamma,gamma,gamma)), fColor.a);
#else
  gl_FragColor = fColor;
#endif
