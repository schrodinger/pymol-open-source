varying vec3 N, L0, H0, L1, H1;
varying vec4 D0, A0, D1, A1;
varying float fog;
void main()
{
  N = normalize(gl_NormalMatrix * gl_Normal);
  vec3 eye_pos = vec3(gl_ModelViewMatrix * gl_Vertex);
  vec3 aux = vec3(gl_LightSource[0].position - vec4(eye_pos, 1.0));
  L0 = normalize(aux);
  H0 = normalize(gl_LightSource[0].halfVector.xyz);
  aux = vec3(gl_LightSource[1].position - vec4(eye_pos, 1.0));
  L1 = normalize(aux);
  H1 = normalize(gl_LightSource[1].halfVector.xyz);
  A0 = gl_LightSource[0].ambient * gl_Color;
  D0 = gl_LightSource[0].diffuse * gl_Color;
  A1 = gl_LightSource[1].ambient * gl_Color;
  D1 = gl_LightSource[1].diffuse * gl_Color;
  gl_FogFragCoord = -eye_pos.z;
  fog = (gl_Fog.end - gl_FogFragCoord) * gl_Fog.scale;
  gl_Position = ftransform();
}
