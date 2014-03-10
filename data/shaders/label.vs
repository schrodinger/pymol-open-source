
attribute vec4 attr_worldpos;
attribute vec3 attr_screenoffset;
attribute vec2 attr_texcoords;
attribute vec3 attr_screenworldoffset;
attribute vec4 attr_t_pickcolor; // changed name for ATI, optional attribute
                                 // cannot be first/0 index, ATI sets index
                                 // alphabetically
varying vec2 textureLookup ;
varying vec3 normalizedViewCoordinate;
varying vec4 pickcolor ;

uniform vec2 screenSize;

uniform float aspectRatioAdjustment;

uniform float screenOriginVertexScale;

varying float fog;
uniform float fog_enabled; // actually bool

void main()
{
  vec4 pos = gl_ModelViewProjectionMatrix * attr_worldpos;
  pos /= pos.w;
  pos.xyz += attr_screenworldoffset / (screenSize.xyx * screenOriginVertexScale);
  pos.xy += attr_screenoffset.xy * 2. / screenSize * vec2(aspectRatioAdjustment, 1.);

  gl_Position = pos;
  textureLookup = attr_texcoords;
  normalizedViewCoordinate = pos.xyz / 2.0 + 0.5;

  if (fog_enabled > 0.5) {
    vec3 eye_pos = vec3(gl_ModelViewMatrix * attr_worldpos);
    fog = max(0.0, (gl_Fog.end - abs(eye_pos.z)) * gl_Fog.scale);
  } else {
    fog = 1.1; // >= 1.0
  }
  pickcolor = attr_t_pickcolor;
}
