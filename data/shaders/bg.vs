attribute vec3 a_Vertex;
varying vec2 bgTextureLookup ;

uniform float isCentered;
uniform vec2 imageSize;
uniform vec2 screenCenter;
uniform vec2 pixelSize;

void main()
{
  vec2 imageVertex = a_Vertex.xy;
  imageVertex = pixelSize * floor(imageVertex / pixelSize);
  gl_Position = vec4(imageVertex.x, imageVertex.y, .5, 1.);
  bgTextureLookup = (1. + a_Vertex.xy) / 2.;
  bgTextureLookup = pixelSize * floor(bgTextureLookup / pixelSize);
}
