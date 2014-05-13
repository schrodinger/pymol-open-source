
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
  pos.xy += attr_screenworldoffset.xy / (screenSize.xy * screenOriginVertexScale);
  pos.xy += attr_screenoffset.xy * 2. / screenSize * vec2(aspectRatioAdjustment, 1.);

  vec3 viewVector = vec3(gl_ModelViewMatrixTranspose * vec4(0.,0., -1.,0.));
  vec4 a_center = (attr_worldpos + attr_screenworldoffset.z * vec4(viewVector, 0.));
  vec4 transformedPositionZ = gl_ModelViewProjectionMatrix * a_center;
  transformedPositionZ.xyz = transformedPositionZ.xyz/transformedPositionZ.w;
  transformedPositionZ.w = 1.;
  pos.z = transformedPositionZ.z;
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
