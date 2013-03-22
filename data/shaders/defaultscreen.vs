attribute vec4 a_Vertex;
attribute vec4 a_Color;

varying vec4 COLOR ;

uniform vec2 t2PixelSize;

void main()
{
  COLOR = a_Color;
  
  gl_Position = vec4(a_Vertex.x * t2PixelSize.x - 1., a_Vertex.y * t2PixelSize.y - 1., .5, 1.);
}
