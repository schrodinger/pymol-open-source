
uniform sampler2D bgTextureMap;
varying vec2 bgTextureLookup ;

uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;

void main()
{
#ifdef bg_image_mode_stretched
  vec2 bgLookup = pixelSize * floor(bgTextureLookup / pixelSize);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  gl_FragColor = vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, 1.);
#endif
#ifdef bg_image_mode_solid
  gl_FragColor = vec4(fogSolidColor.rgb, 1.);
#endif
}

