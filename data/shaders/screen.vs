#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

attribute vec3 attr_screenoffset;
attribute vec2 attr_texcoords;
attribute vec4 attr_backgroundcolor;

varying vec2 textureLookup ;
varying vec4 backgroundColor;
uniform vec2 t2PixelSize;

void main()
{
  vec4 transformedPosition = vec4(attr_screenoffset.x * t2PixelSize.x - 1.,
                                  attr_screenoffset.y * t2PixelSize.y - 1.,
                                  .9, 1.);
  gl_Position = transformedPosition;
  backgroundColor = attr_backgroundcolor;
  textureLookup = attr_texcoords;
}
