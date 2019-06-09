#include webgl_header.fs

#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

#include anaglyph_header.fs

varying vec2 bgTextureLookup ;

#include compute_fog_color.fs

void main()
{
  gl_FragColor = vec4(ComputeBgColor(), 1.);
}

