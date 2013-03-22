
varying vec4 COLOR ;

#include ANAGLYPH_HEADER

void main()
{
  vec4 fColor = COLOR;

#include ANAGLYPH_BODY
}

