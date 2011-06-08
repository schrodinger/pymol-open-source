varying vec3 N, L0, H0, L1, H1;
varying vec4 D0, A0, D1, A1;

uniform float fog_r;
uniform float fog_g;
uniform float fog_b;
uniform float fog_enabled;

varying float fog;

void main()
{
  vec3 n, h;
  float NdotL, NdotH;
  vec4 color = A0 + A1;
  n = normalize(N);
  NdotL = max(dot(n, normalize(L0)), 0.0);
  float shininess = gl_FrontMaterial.shininess;
  if (NdotL > 0.0) {
      color += D0 * NdotL;
      h = normalize(H0);
      NdotH = max(dot(n, h), 0.0);
      color += gl_LightSource[0].specular * pow(NdotH, shininess);
  }
  NdotL = max(dot(n, normalize(L1)), 0.0);
  if (NdotL > 0.0) {
      color += D1 * NdotL;
      h = normalize(H1);
      NdotH = max(dot(n, h), 0.0);
      color += gl_LightSource[1].specular * pow(NdotH, shininess);
  }

  vec3 fog_color = vec3(fog_r, fog_g, fog_b);
  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);

  gl_FragColor = vec4(mix(fog_color, color.rgb, cfog), color.a);
}

