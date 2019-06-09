#include webgl_header.fs

#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

/* ios recent */
uniform sampler2D textureMap;
varying vec2 textureLookup ;
varying vec4 pickcolor ;

varying float fog;
varying vec3 normalizedViewCoordinate;
#define bgTextureLookup normalizedViewCoordinate.xy

#include anaglyph_header.fs
#include compute_fog_color.fs

void main()
{
  if (isPicking) {
    gl_FragColor = pickcolor;
  } else {
    gl_FragColor = texture2D(textureMap, textureLookup);

#ifdef PYMOL_IOS
    if (gl_FragColor.a < .1)
      discard;
#else
    if (gl_FragColor.a < .05)
      discard;
#endif

    gl_FragColor = ApplyColorEffects(gl_FragColor, gl_FragCoord.z);
    gl_FragColor = ApplyFog(gl_FragColor, fog);

#ifdef NO_ORDER_TRANSP
    setFragDataForNoOrderTransp(gl_FragCoord.z);
#endif
  }
}

