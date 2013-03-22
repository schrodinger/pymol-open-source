attribute vec3 attr_screenoffset;
attribute vec2 attr_texcoords;
varying vec2 textureLookup ;

uniform vec2 t2PixelSize;


void main()
{
  vec4 transformedPosition = vec4(attr_screenoffset.x * t2PixelSize.x - 1.,
                                  attr_screenoffset.y * t2PixelSize.y - 1.,
                                  .9, 1.);
  gl_Position = transformedPosition;
  textureLookup = attr_texcoords;
}
