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
struct FogParameters {
  vec4 color;
  float density;
  float start;
  float end;
  float scale;
};
uniform FogParameters st_Fog;
uniform mat4 u_ModelViewMatrix;
uniform mat3 u_NormalMatrix;
attribute vec4 a_Vertex;
attribute vec3 a_Normal;
attribute vec4 a_Color;
varying vec3 N, L0, H0, L1, H1;
varying vec4 D0, A0, D1, A1;
varying float fog;
varying float v_FogFragCoord;
void main()
{
  N = normalize(u_NormalMatrix * a_Normal);
  vec3 eye_pos = vec3(u_ModelViewMatrix * a_Vertex);
  vec3 aux = vec3(u_LightSource[0].position - vec4(eye_pos, 1.0));
  L0 = normalize(aux);
  H0 = normalize(u_LightSource[0].halfVector.xyz);
  aux = vec3(u_LightSource[1].position - vec4(eye_pos, 1.0));
  L1 = normalize(aux);
  H1 = normalize(u_LightSource[1].halfVector.xyz);
  A0 = u_LightSource[0].ambient * a_Color;
  D0 = u_LightSource[0].diffuse * a_Color;
  A1 = u_LightSource[1].ambient * a_Color;
  D1 = u_LightSource[1].diffuse * a_Color;
  v_FogFragCoord = -eye_pos.z;
  fog = (st_Fog.end - v_FogFragCoord) * st_Fog.scale;
  gl_Position =  u_ModelViewMatrix * a_Vertex; // this was converted from ftransform(), is that right?
}
