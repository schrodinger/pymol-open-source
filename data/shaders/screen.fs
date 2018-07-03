#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

/* ios recent */
uniform sampler2D textureMap;

varying vec4 backgroundColor;
varying vec2 textureLookup ;

void main()
{
  vec4 fColor = texture2D(textureMap, textureLookup);
  vec4 bColor = backgroundColor * (1.-fColor.a);
  gl_FragColor = vec4(bColor.rgb + fColor.rgb * fColor.a, bColor.a + fColor.a );

#ifdef ANAGLYPH
  gl_FragColor.rgb = vec3((gl_FragColor.r + gl_FragColor.b + gl_FragColor.g) / 3.);
#endif
}

