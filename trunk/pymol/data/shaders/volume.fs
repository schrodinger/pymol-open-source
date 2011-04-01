uniform sampler3D volumeTex;
uniform sampler1D colorTex;
varying float fog;

void main()
{
  vec3 fog_color = vec3(0.0, 0.0, 0.0);
  vec4 color = texture1D(colorTex, texture3D(volumeTex, gl_TexCoord[0].xyz).x);
  if (color.a == 0.0) discard;
  gl_FragColor = vec4(vec3(mix(fog_color, color.rgb, fog)), color.a);
}

