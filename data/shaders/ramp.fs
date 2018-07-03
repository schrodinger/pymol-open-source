#ifdef PURE_OPENGL_ES_2
precision lowp float;
#endif
varying vec4 COLOR ;
varying vec3 NORMAL ;

void main()
{
  float NdotV = dot(NORMAL, vec3(0.,0.,1.));
  gl_FragColor = vec4(NdotV * COLOR.xyz, COLOR.a);

#ifdef ANAGLYPH
  gl_FragColor.rgb = vec3((gl_FragColor.r + gl_FragColor.b + gl_FragColor.g) / 3.);
#endif
}

