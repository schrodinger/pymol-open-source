#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

#include webgl_header.vs

attribute vec4 attr_targetpos;
attribute vec4 attr_worldpos;
attribute vec3 attr_screenoffset;
attribute vec2 attr_texcoords;
attribute vec3 attr_screenworldoffset;
attribute vec4 attr_pickcolor;
attribute float attr_relative_mode;

varying vec2 textureLookup ;
varying vec3 normalizedViewCoordinate;
varying vec4 pickcolor ;

uniform vec2 screenSize;

uniform float screenOriginVertexScale;

uniform float scaleByVertexScale;
uniform float labelTextureSize;
varying float fog;

uniform float front;
uniform float clipRange;

vec4 normalizeVec4(vec4 point){
     return vec4(point.xyz/point.w, 1.);
}

float convertNormalZToScreenZ(float normalz){
   float a_centerN = (normalz + 1.) / 2.;
   float ptInPreProjectionZ = -(front + clipRange * a_centerN);
   vec4 ptInPreProjection = vec4(0., 0., ptInPreProjectionZ, 1.);
   vec4 projVect = g_ProjectionMatrix * ptInPreProjection;
   return projVect.z / projVect.w;
}

void main()
{
  float drawConnector, isProjected, isScreenCoord, isPixelCoord, zTarget;

  drawConnector = step(1., mod(attr_relative_mode, 2.));
  isScreenCoord = step(2., mod(attr_relative_mode, 4.));
  isPixelCoord = step(4., mod(attr_relative_mode, 8.));
  zTarget = step(8., mod(attr_relative_mode, 16.));
  isProjected = step(isPixelCoord + isScreenCoord, 0.5);

  vec3 viewVector = vec3(vec4(0.,0.,-1.,0.) * g_ModelViewMatrix);
  float sovx = screenOriginVertexScale;

#ifdef openvr_enabled
  // calc real view vector for label (not negative z-axis in camera space)
  vec3 camPos = vec4(gl_ModelViewMatrixInverse * vec4(0.,0.,0.,1.)).xyz;
  viewVector = normalize(attr_worldpos.xyz - camPos);
  // don't use vertexscale for openvr, because the scale is calced from Pos.z as a distance 
  sovx = 1.0 / 2.0;
#endif

  float screenVertexScale = scaleByVertexScale * sovx * labelTextureSize + (1. - scaleByVertexScale);

  vec4 transformedPosition = g_ProjectionMatrix * g_ModelViewMatrix * attr_worldpos;
  vec4 targetPosition = normalizeVec4(g_ProjectionMatrix * g_ModelViewMatrix * attr_targetpos);
  transformedPosition.xyz = transformedPosition.xyz/transformedPosition.w;
  transformedPosition.xy = ( floor(transformedPosition.xy * screenSize + .5 ) + .5 ) / screenSize;

  vec4 a_center = (attr_worldpos + attr_screenworldoffset.z * vec4(viewVector, 0.));
  vec4 transformedPositionZ = g_ProjectionMatrix * g_ModelViewMatrix * a_center;
  transformedPositionZ.xyz = transformedPositionZ.xyz/transformedPositionZ.w;
  transformedPositionZ.w = 1.;
  vec2 pixOffset = ((2. * attr_worldpos.xy / screenSize) - 1.);
  transformedPosition = isProjected * transformedPosition + isScreenCoord * attr_worldpos + isPixelCoord * vec4(pixOffset.x, pixOffset.y, -0.5, 0.);
  transformedPosition.xy = transformedPosition.xy + attr_screenworldoffset.xy/(screenSize*sovx);

  transformedPosition.z = (1.-zTarget) * ((isProjected * transformedPositionZ.z) + (1.-isProjected) * convertNormalZToScreenZ(attr_worldpos.z)) + zTarget * targetPosition.z;

  transformedPosition.xy += attr_screenoffset.xy * 2./(screenSize.xy*screenVertexScale);

  transformedPosition.w = 1.;
  gl_Position = transformedPosition;
  textureLookup = attr_texcoords;
  normalizedViewCoordinate = (gl_Position.xyz/gl_Position.w) / 2.0 + 0.5;

#ifdef depth_cue
  vec3 eye_pos = mix(attr_worldpos, g_ModelViewMatrix * attr_worldpos, isProjected).xyz;
  fog = (g_Fog_end + eye_pos.z) * g_Fog_scale;
    fog = max(fog, 0.0);
#else
    fog = 1.1; // >= 1.0
#endif

  pickcolor = attr_pickcolor;
}
