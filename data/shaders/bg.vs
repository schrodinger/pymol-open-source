#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

attribute vec3 a_Vertex;
varying vec2 bgTextureLookup ;

uniform bool isPicking;

void main()
{
  gl_Position = vec4(a_Vertex.xy, .5, 1.);
  bgTextureLookup = (1. + a_Vertex.xy) / 2.;
}
