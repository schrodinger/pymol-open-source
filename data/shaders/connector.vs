
attribute vec4 a_target_pt3d;
attribute vec4 a_center_pt3d;
attribute vec3 a_indentFactor;
attribute vec3 a_screenWorldOffset;
attribute vec2 a_textSize;
attribute vec4 a_Color;
attribute float a_relative_mode;
attribute float a_draw_flags;
attribute vec4 a_bkgrd_color;
attribute float a_rel_ext_length;
attribute float a_con_width;

varying vec4 COLORIn;

#ifdef use_geometry_shaders
varying vec4 a_center;
varying vec4 a_target;
varying vec2 a_textSizeIn;
varying vec3 a_screenWorldOffsetIn;
varying vec2 a_indentFactorIn;
varying float a_relative_modeIn;
varying float a_draw_flagsIn;
varying vec4 a_bkgrd_colorIn;
varying float a_rel_ext_lengthIn;
varying float a_con_widthIn;
varying vec4 COLOR;
varying float lineEdge;
varying float aaCutoff;
#else
attribute float a_isCenterPt;
#endif

varying vec3 NORMAL;
varying float fog;
varying vec2 bgTextureLookup;

#define E 2.718281828459045

uniform vec2 screenSize;
uniform float screenOriginVertexScale;
uniform float labelTextureSize;
uniform float front;
uniform float clipRange;
float PI = 3.14159265358979323846264;

#include connector.shared


vec4 getGLPosition(vec4 endpointOnBBX, vec4 transformedTarget, vec4 endpointExtendedOffBBX, float connector_mode_2_4){
#ifdef use_geometry_shaders
  float a_isCenterPt = 1.;
#endif
  float isExtended = step(1.5, a_isCenterPt);
  float isTarget = (1.-isExtended) * step(.5, a_isCenterPt );
  float isCenter = (1.-isTarget) * (1.-isExtended);
  return connector_mode_2_4 * ( (isCenter * endpointOnBBX) + (isTarget * transformedTarget) + (isExtended * endpointExtendedOffBBX)) + 
    	 (1.-connector_mode_2_4) * ( (step(.5,a_isCenterPt) * endpointOnBBX) + (1.-step(.5,a_isCenterPt)) * transformedTarget);
}

void main()
{
#ifdef use_geometry_shaders
  a_center = a_center_pt3d ;
  a_target = a_target_pt3d ;
  a_textSizeIn = a_textSize;
  a_screenWorldOffsetIn = a_screenWorldOffset;
  a_indentFactorIn = a_indentFactor.xy;
  gl_Position = a_center_pt3d ;
  a_relative_modeIn = a_relative_mode;
  a_draw_flagsIn = a_draw_flags;
  a_bkgrd_colorIn = a_bkgrd_color;
  a_rel_ext_lengthIn = a_rel_ext_length;
  a_con_widthIn = a_con_width;
#else
  NORMAL = vec3(0.,0.,1.);

  vec4 transformedCenter, transformedTarget;
  vec2 textSizeInScreen, offset, drawVector;
  float doNotDraw, isProjected, zValue, zTarget;
  float connector_mode_0, connector_mode_1, connector_mode_2, connector_mode_3, connector_mode_4, drawBackgroundOutline, drawBackground, drawConnector, connector_mode_2_4 ;
  getDrawFlags(a_draw_flags, connector_mode_0, connector_mode_1, connector_mode_2, connector_mode_3, connector_mode_4, drawBackgroundOutline, drawBackground, drawConnector);
  connector_mode_2_4 = connector_mode_2 + connector_mode_4;
  calculatePreConnectorInfo(a_center_pt3d, a_target_pt3d, a_textSize, a_indentFactor.xy, a_screenWorldOffset, transformedCenter, transformedTarget, textSizeInScreen, offset, drawVector, doNotDraw, a_relative_mode, isProjected, zValue, a_con_width, zTarget);
  vec4 endpointOnBBX, endpointExtendedOffBBX;
  calculateConnectorInfo(a_center_pt3d, textSizeInScreen, a_textSize, drawVector, offset, fog, transformedCenter, transformedTarget, endpointOnBBX, endpointExtendedOffBBX, connector_mode_0, connector_mode_1, connector_mode_2, connector_mode_3, connector_mode_4, a_rel_ext_length, isProjected, zValue);

  // if the zValue is not within the clipping planes, then set the transformedTarget z
  // to the zValue as well, to make sure that the connector doesn't get drawn
  float withinView = step(zValue, 1.) * step(-1., zValue);
  transformedTarget.z = withinView * transformedTarget.z + (1.-withinView) * zValue;

  gl_Position = (1.-doNotDraw) * getGLPosition(endpointOnBBX, transformedTarget, endpointExtendedOffBBX, connector_mode_2_4);
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
#endif

  COLORIn = a_Color;

}
