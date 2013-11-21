attribute vec4 attr_worldpos;
attribute vec3 attr_screenoffset;
attribute vec2 attr_texcoords;
attribute vec3 attr_screenworldoffset;
attribute vec4 attr_pickcolor;
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
  vec4 transformedPosition = gl_ModelViewProjectionMatrix * attr_worldpos;
  transformedPosition.xyz = transformedPosition.xyz/transformedPosition.w;
  transformedPosition.w = 1.;
  transformedPosition.xy = transformedPosition.xy + attr_screenworldoffset.xy/(screenSize*screenOriginVertexScale);
  transformedPosition.z = transformedPosition.z + attr_screenworldoffset.z/(screenSize.x*screenOriginVertexScale);
  transformedPosition.x = transformedPosition.x + aspectRatioAdjustment * attr_screenoffset.x * 2./screenSize.x;
  transformedPosition.y = transformedPosition.y + attr_screenoffset.y * 2./screenSize.y;
  gl_Position = transformedPosition;
  textureLookup = attr_texcoords;
  normalizedViewCoordinate = (gl_Position.xyz/gl_Position.w) / 2.0 + 0.5;
  if (fog_enabled > 0.5) {
    vec3 eye_pos = vec3(gl_ModelViewMatrix * attr_worldpos);
    fog = max(0.0, (gl_Fog.end - abs(eye_pos.z)) * gl_Fog.scale);
  } else {
    fog = 1.1; // >= 1.0
  }
  pickcolor = attr_pickcolor;
}
