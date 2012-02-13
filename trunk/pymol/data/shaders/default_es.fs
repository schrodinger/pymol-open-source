precision mediump float;
const int a_MaxLights=4;
struct st_LightSourceParameters { 
  vec4 ambient;
  vec4 diffuse;
  vec4 specular;
  vec4 position;
  vec4 halfVector;
  vec3 spotDirection;
  float spotExponent;
  float spotCutoff;
  float spotCosCutoff;
  float constantAttenuation;
  float linearAttenuation;
  float quadraticAttenuation;
}; 
uniform st_LightSourceParameters u_LightSource[a_MaxLights];
struct gl_MaterialParameters {
  vec4 emission;
  vec4 ambient;
  vec4 diffuse;
  vec4 specular;
  float shininess;
};
uniform gl_MaterialParameters st_FrontMaterial;
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
  float shininess = st_FrontMaterial.shininess;
  if (NdotL > 0.0) {
      color += D0 * NdotL;
      h = normalize(H0);
      NdotH = max(dot(n, h), 0.0);
      color += u_LightSource[0].specular * pow(NdotH, shininess);
  }
  NdotL = max(dot(n, normalize(L1)), 0.0);
  if (NdotL > 0.0) {
      color += D1 * NdotL;
      h = normalize(H1);
      NdotH = max(dot(n, h), 0.0);
      color += u_LightSource[1].specular * pow(NdotH, shininess);
  }

  vec3 fog_color = vec3(fog_r, fog_g, fog_b);
  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);

  gl_FragColor = vec4(mix(fog_color, color.rgb, cfog), color.a);
}
