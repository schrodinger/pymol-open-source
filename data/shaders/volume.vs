#include webgl_header.vs

uniform vec3 slabTexOffset;

varying vec3 vertexM;

varying float fog;
varying vec2 bgTextureLookup;

void main()
{
  vec4 vertex = g_ModelViewMatrix * gl_Vertex;
  gl_TexCoord[0] = gl_MultiTexCoord0;

#ifdef volume_mode
#endif // volume_mode

  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * gl_Vertex;
  fog = (g_Fog_end + vertex.z) * g_Fog_scale;
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
