uniform sampler3D volumeTex;
uniform sampler1D colorTex;

uniform float fog_r;
uniform float fog_g;
uniform float fog_b;
uniform float fog_enabled;
varying float fog;

void main()
{
  vec3 fog_color = vec3(fog_r, fog_g, fog_b);
  vec4 color = texture1D(colorTex, texture3D(volumeTex, gl_TexCoord[0].xyz).x);
  if (color.a == 0.0) discard;
  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  gl_FragColor = vec4(vec3(mix(fog_color, color.rgb, cfog)), color.a);
}

