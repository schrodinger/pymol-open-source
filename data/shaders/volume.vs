varying float fog;
uniform float g_Fog_end;
uniform float g_Fog_scale;
varying vec2 bgTextureLookup;

void main()
{
  vec4 vertex = gl_ModelViewMatrix * gl_Vertex;
  gl_TexCoord[0] = gl_MultiTexCoord0;
  gl_ClipVertex = vertex;
  gl_Position = ftransform();
  gl_FogFragCoord = -vertex.z;
  fog = (g_Fog_end - gl_FogFragCoord) * g_Fog_scale;
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
