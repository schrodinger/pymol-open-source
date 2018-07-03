#version 120
#extension GL_EXT_geometry_shader4 : require
#extension GL_EXT_gpu_shader4 : require

varying out vec3 NORMAL ;
varying out vec4 COLOR ;
varying out float fog;
varying out vec2 bgTextureLookup;
varying out float lineEdge;
varying out float aaCutoff;

varying in vec4 a_center[];
varying in vec4 a_target[];
varying in vec4 COLORIn[];
varying in vec2 a_textSizeIn[];
varying in vec3 a_screenWorldOffsetIn[];
varying in vec2 a_indentFactorIn[];
varying in float a_relative_modeIn[];
varying in float a_draw_flagsIn[];
varying in vec4 a_bkgrd_colorIn[];
varying in float a_rel_ext_lengthIn[];
varying in float a_con_widthIn[];

uniform vec2 screenSize;
uniform float screenOriginVertexScale;
uniform float antialiasedLines;
uniform float textureToLabelSize;
uniform float front;
uniform float clipRange;

float PI = 3.14159265358979323846264;

float mysmoothstep(float edge0, float edge1, float x){
  float rets = step(edge0, edge1) * step(edge1, edge0);
  return rets * step(edge0, x) + (1.-rets) * smoothstep(edge0, edge1, x);
}

#include connector.shared

struct lineSeg_t {
  vec2 p1;
  vec2 p2;
};
bool lineSegIntersection(vec2 p1,  vec2 p2, vec2 p3, vec2 p4, out vec2 po) {
	float  distAB, theCos, theSin, newX, abpos;
	vec2 p2i, p3i, p4i;
	if ((p1.x==p2.x && p1.y==p2.y) || (p3.x==p4.x && p3.y==p4.y)) return false;
	if ((p1.x==p3.x && p1.y==p3.y) || (p2.x==p3.x && p2.y==p3.y)
	||  (p1.x==p4.x && p1.y==p4.y) || (p2.x==p4.x && p2.y==p4.y)) {
	return false; }
	p2i.x = p2.x; p2i.y = p2.y;
	p3i.x = p3.x; p3i.y = p3.y;
	p4i.x = p4.x; p4i.y = p4.y;
	p2i.x -= p1.x; p2i.y -= p1.y;
	p3i.x -= p1.x; p3i.y -= p1.y;
	p4i.x -= p1.x; p4i.y -= p1.y;
	distAB = sqrt(p2i.x * p2i.x + p2i.y * p2i.y);
	theCos = p2i.x / distAB;
	theSin = p2i.y / distAB;
	newX = p3i.x * theCos + p3i.y * theSin;
	p3i.y = p3i.y * theCos - p3i.x * theSin;
	p3i.x = newX;
	newX = p4i.x * theCos + p4i.y * theSin;
	p4i.y = p4i.y * theCos - p4i.x * theSin;
	p4i.x = newX;
	if ((p3i.y < 0. && p4i.y<0.) || (p3i.y >= 0. && p4i.y >= 0.)) return false;
	abpos = p4i.x + (p3i.x - p4i.x) * p4i.y / (p4i.y - p3i.y);
	if (abpos<0. || abpos>distAB) return false;
	po.x = p1.x + abpos * theCos;
	po.y = p1.y + abpos * theSin;
	return true;
}

void setBG(){
	bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}

void drawLineAsGeometryClipped(vec4 pt1, vec4 pt2, lineSeg_t labelTop, lineSeg_t labelBottom, lineSeg_t labelLeft, lineSeg_t labelRight, out vec2 dirv, float lineWidth, float aa){
	vec4 pt1E = pt1;
	vec4 pt2E = pt2;
	pt1E.xy = floor(pt1.xy * screenSize) / screenSize;
	pt2E.xy = floor(pt2.xy * screenSize) / screenSize;
	vec3 diffV = pt1E.xyz - pt2E.xyz;
	vec2 dirvPnorm = vec2(normalize(cross(diffV, vec3(0.,0.,1.))).xy);
	vec2 dirvP = lineWidth*dirvPnorm;
	dirv = dirvP/screenSize;
        int l1V = 0 , l2V = 0;
        lineSeg_t l1, l2;
        vec2 isec1, isec2;
        l1.p1 = pt1.xy + dirv;
        l1.p2 = pt2.xy + dirv;
        l2.p1 = pt1.xy - dirv;
        l2.p2 = pt2.xy - dirv;

        if (lineSegIntersection(labelTop.p1, labelTop.p2, l1.p1, l1.p2, isec1)) {
              l1V = 1;
        } else if (lineSegIntersection(labelRight.p1, labelRight.p2, l1.p1, l1.p2, isec1)) {
              l1V = 2;
        } else if (lineSegIntersection(labelBottom.p1, labelBottom.p2, l1.p1, l1.p2, isec1)) {
              l1V = 3;
        } else if (lineSegIntersection(labelLeft.p1, labelLeft.p2, l1.p1, l1.p2, isec1)) {
              l1V = 4;
        } 

        if (lineSegIntersection(labelTop.p1, labelTop.p2, l2.p1, l2.p2, isec2)) {
              l2V = 1;
        } else if (lineSegIntersection(labelRight.p1, labelRight.p2, l2.p1, l2.p2, isec2)) {
              l2V = 2;
        } else if (lineSegIntersection(labelBottom.p1, labelBottom.p2, l2.p1, l2.p2, isec2)) {
              l2V = 3;
        } else if (lineSegIntersection(labelLeft.p1, labelLeft.p2, l2.p1, l2.p2, isec2)) {
              l2V = 4;
        } 

        // Both lines run through one side
        if ( (l1V > 0 && l2V > 0) && (l1V == l2V) ) {
              pt1E.xy = pt1.xy + dirv; pt1E.z = pt1.z;
              gl_Position = pt1E;
              lineEdge = 1.; aaCutoff = aa;
              setBG(); EmitVertex();

              pt1E.xy = pt1.xy - dirv; pt1E.z = pt1.z;
              gl_Position = pt1E;
              lineEdge = -1.;
              setBG(); EmitVertex();
              
              pt2E.xy = isec1;  pt2E.z = pt2.z;
              gl_Position = pt2E;
              lineEdge = 1.;
              setBG(); EmitVertex();

              pt2E.xy = isec2;  pt2E.z = pt2.z;
              gl_Position = pt2E;
              lineEdge = -1.;
              setBG(); EmitVertex();
              EndPrimitive();
              lineEdge = 0.;
              aaCutoff = 0.; 

        } else if ( (l1V > 0 && l2V > 0) && (l1V != l2V) ) { // Find the corner 
              vec2 isecc;
              if ( (l1V == 1 && l2V == 2) || (l1V == 2 && l2V == 1) ) { // UR
                      isecc = labelTop.p1;
              } else if ( (l1V == 2 && l2V == 3) || (l1V == 3 && l2V == 2) ) { // LR
                      isecc = labelBottom.p2;
              } else if ( (l1V == 3 && l2V == 4) || (l1V == 4 && l2V == 3) ) { // LL
                      isecc = labelBottom.p1;
              } else if ( (l1V == 4 && l2V == 1) || (l1V == 1 && l2V == 4) ) { // UL
                      isecc = labelTop.p2;
              }
              
              pt2E.xy = isec1.xy; pt2E.z = pt2.z;
              gl_Position = pt2E;
              lineEdge = 1.; aaCutoff = aa;
              setBG(); EmitVertex();

              pt1E.xy = pt1.xy + dirv; pt1E.z = pt1.z;
              gl_Position = pt1E;
              lineEdge = 1.;
              setBG(); EmitVertex();

              pt2E.xy = isecc.xy; pt2E.z = pt2.z;
              gl_Position = pt2E;


              /* since we are adding a corner vertex, we also need to 
                 account for anti-aliasing by setting lineEdge for that 
                 vertex. */              
              float len1 = length(screenSize* (isecc.xy - isec1.xy)); // length between corner and first intersection
              float len2 = length(screenSize* (isecc.xy - isec2.xy)); // length between corner and second intersection
              float lenc = length(screenSize* (isec2.xy - isec1.xy)); // length between intersections, i.e., line width
              float ang = atan(len2 / len1);  // angle between perpendicular vector and first edge
              float n2 = sin(ang) * len2;     // second edge projected onto perpendicular vector
              float lenr = (n2 / lenc) * 2. - 1.;  // figure out ratio and normalize it between -1 and 1
              float lenrsign = 2. * (step(0., lenr) - .5);  // sign for lenr used below
              lineEdge = lenrsign * (1. - mysmoothstep(0., aa, 1 - abs(lenr))); // take into account aa for corner vertex lineEdge

              setBG(); EmitVertex();

              pt1E.xy = pt1.xy - dirv; pt1E.z = pt1.z;
              gl_Position = pt1E;
              lineEdge = -1.;
              setBG(); EmitVertex();

              pt2E.xy = isec2.xy; pt2E.z = pt2.z;
              gl_Position = pt2E;
              lineEdge = -1.;
              setBG(); EmitVertex();
              EndPrimitive();
              
         
        } 
}

void drawLineAsGeometry(vec4 pt1, vec4 pt2, out vec2 dirv, float lineWidth, float aa){
	vec4 pt1E = pt1;
	vec4 pt2E = pt2;
	pt1E.xy = floor(pt1.xy * screenSize) / screenSize;
	pt2E.xy = floor(pt2.xy * screenSize) / screenSize;
	vec3 diffV = pt1E.xyz - pt2E.xyz;
	vec2 dirvPnorm = vec2(normalize(cross(diffV, vec3(0.,0.,1.))).xy);
	vec2 dirvP = lineWidth*dirvPnorm;
	dirv = dirvP/screenSize;

	pt1E.xy = pt1.xy + dirv;
	gl_Position = pt1E;
	lineEdge = -1.; aaCutoff = aa;
	setBG(); EmitVertex();

	pt2E.xy = pt2.xy + dirv;
	gl_Position = pt2E;
	lineEdge = -1.;
	setBG(); EmitVertex();

	pt1E.xy = pt1.xy - dirv;
	gl_Position = pt1E;
	lineEdge = 1.;
	setBG(); EmitVertex();

	pt2E.xy = pt2.xy - dirv;
	gl_Position = pt2E;
	lineEdge = 1.;
	setBG(); EmitVertex();
	EndPrimitive();
	lineEdge = 0.;
	aaCutoff = 0.; 
}

void drawLineAsGeometryWithOffsets(vec4 pt1, vec4 pt2, float topext, float bottomext){
	vec4 pt1E = pt1;
	vec4 pt2E = pt2;
	vec3 diffV = pt1E.xyz - pt2E.xyz;
	vec2 linev = vec2(a_con_widthIn[0] * normalize(diffV).xy) / screenSize;
	vec2 dirv = vec2((a_con_widthIn[0]*normalize(cross(diffV, vec3(0.,0.,1.)).xy)))/screenSize;
	pt1E.xy = pt1.xy + dirv + topext * linev;
	gl_Position = pt1E;
	setBG(); EmitVertex();

	pt2E.xy = pt2.xy + dirv - topext * linev;
	gl_Position = pt2E;
	setBG(); EmitVertex();

	pt1E.xy = pt1.xy - dirv + bottomext * linev;
	gl_Position = pt1E;
	setBG(); EmitVertex();

	pt2E.xy = pt2.xy - dirv - bottomext * linev;
	gl_Position = pt2E;
	setBG(); EmitVertex();
	EndPrimitive();
}

void main()
{
	vec4 transformedCenter, transformedTarget;
	vec2 textSizeInScreen, offset, drawVector;
	float doNotDraw;
	float connector_mode_0, connector_mode_1, connector_mode_2, connector_mode_3, connector_mode_4, connector_mode_2_4, drawBackgroundOutline, drawBackground, drawConnector, isProjected, zValue, zTarget;
	lineEdge = 0.;
	aaCutoff = 0.;
	getDrawFlags(a_draw_flagsIn[0], connector_mode_0, connector_mode_1, connector_mode_2, connector_mode_3, connector_mode_4, drawBackgroundOutline, drawBackground, drawConnector);

	float negLabelSize = step(textureToLabelSize, 0.);
	float textureToLabelSizeA = abs(textureToLabelSize);
	float extLength = (1.-negLabelSize) * a_rel_ext_lengthIn[0] + negLabelSize * (-a_rel_ext_lengthIn[0]/(screenOriginVertexScale*2.) );
	vec2 textSize = a_textSizeIn[0] / textureToLabelSizeA;
	connector_mode_2_4 = connector_mode_2 + connector_mode_4;
	calculatePreConnectorInfo(a_center[0], a_target[0], textSize, a_indentFactorIn[0], a_screenWorldOffsetIn[0], transformedCenter, transformedTarget, textSizeInScreen, offset, drawVector, doNotDraw, a_relative_modeIn[0], isProjected, zValue, a_con_widthIn[0], zTarget);
	vec4 endpointOnBBX, endpointExtendedOffBBX;
	calculateConnectorInfo(a_center[0], textSizeInScreen, textSize, drawVector, offset, fog, transformedCenter, transformedTarget, endpointOnBBX, endpointExtendedOffBBX, connector_mode_0, connector_mode_1, connector_mode_2, connector_mode_3, connector_mode_4, extLength, isProjected, zValue);
	NORMAL = vec3(0.,0.,1.);

        // if the zValue is not within the clipping planes, then set the transformedTarget z
        // to the zValue as well, to make sure that the connector doesn't get drawn
        float withinView = step(zValue, 1.) * step(-1., zValue);
        float ttz = (1.-zTarget) * (1.-isProjected); // not zTarget and not projected, then keep same, otherwise zValue
        transformedTarget.z =  ttz * transformedTarget.z + (1.-ttz) * zValue;

	vec4 transformedPosition3 = transformedCenter;
	transformedPosition3.xy += offset;
	transformedPosition3.z = zValue;
	// bbx of the label
	vec4 upperLeft = transformedPosition3;
	upperLeft.xy += textSizeInScreen * vec2(-1.,1.);
	vec4 upperRight = transformedPosition3;
	upperRight.xy += textSizeInScreen * vec2(1.,1.);
	vec4 lowerLeft = transformedPosition3;
	lowerLeft.xy +=  textSizeInScreen * vec2(-1.,-1.);
	vec4 lowerRight = transformedPosition3;
	lowerRight.xy += textSizeInScreen * vec2(1.,-1.);

        lineSeg_t labelTop, labelBottom, labelLeft, labelRight;

        labelTop.p1 = upperRight.xy;
        labelTop.p2 = upperLeft.xy;
        labelBottom.p1 = lowerLeft.xy;
        labelBottom.p2 = lowerRight.xy;
        labelRight.p1 = upperRight.xy;
        labelRight.p2 = lowerRight.xy;
        labelLeft.p1 = upperLeft.xy;
        labelLeft.p2 = lowerLeft.xy;

	if (drawBackgroundOutline + drawBackground > 0.){
		if (drawBackground > 0.){
                   if (drawBackgroundOutline > 0.){
			vec2 lineWidth = a_con_widthIn[0] / (screenSize);
			vec4 upperLeftM = upperLeft + vec4(lineWidth.x, -lineWidth.y, 0., 0.),
			upperRightM = upperRight + vec4(-lineWidth.x, -lineWidth.y, 0., 0.),
			lowerLeftM = lowerLeft + vec4(lineWidth.x, lineWidth.y, 0., 0.),
			lowerRightM = lowerRight + vec4(-lineWidth.x, lineWidth.y, 0., 0.);

			COLOR = a_bkgrd_colorIn[0];
			gl_Position = upperLeftM; setBG(); EmitVertex();
			gl_Position = upperRightM; setBG(); EmitVertex();
			gl_Position = lowerLeftM;  setBG(); EmitVertex();
			gl_Position = lowerRightM; setBG(); EmitVertex();
			EndPrimitive();
                   } else {
			COLOR = a_bkgrd_colorIn[0];
			gl_Position = upperLeft; setBG(); EmitVertex();
			gl_Position = upperRight; setBG(); EmitVertex();
			gl_Position = lowerLeft;  setBG(); EmitVertex();
			gl_Position = lowerRight; setBG(); EmitVertex();
			EndPrimitive();
                   }
		}
		if (drawBackgroundOutline > 0.){
			COLOR = COLORIn[0];
			drawLineAsGeometryWithOffsets(upperLeft, upperRight, 1., 0.);
			drawLineAsGeometryWithOffsets(lowerRight, upperRight, 0., 1.);
			drawLineAsGeometryWithOffsets(lowerLeft, lowerRight, 0., 1.);
			drawLineAsGeometryWithOffsets(lowerLeft, upperLeft, 1., 0.);
		}
	}
	if (drawConnector > 0 && doNotDraw == 0.){
		COLOR = COLORIn[0];
		float aa = min(1., 2./a_con_widthIn[0]);
		float widthAdjustment = aa * a_con_widthIn[0];
		float adjLineWidth = antialiasedLines * widthAdjustment + a_con_widthIn[0];
		if (connector_mode_2_4 > 0.){
			// create 2 lines, between transformedTarget -> endpointExtendedOffBBX, and endpointExtendedOffBBX -> endpointOnBBX
			vec2 dirv1, dirv2;
			drawLineAsGeometry(endpointExtendedOffBBX, endpointOnBBX, dirv2, a_con_widthIn[0], 0);
			drawLineAsGeometry(transformedTarget, endpointExtendedOffBBX, dirv1, adjLineWidth, aa);
			float norder = step(0., cross(vec3(dirv1.x, dirv1.y, 0.), vec3(dirv2.x, dirv2.y, 0.)).z);
			// need to draw 'fan' that creates a curved elbow between the lines

			lineEdge = 0.; aaCutoff = 0.;
			gl_Position = endpointExtendedOffBBX;
			setBG(); EmitVertex();
			lineEdge = -1.; aaCutoff = aa;
			if (norder < .5){
				vec4 pt = endpointExtendedOffBBX;
				pt.xy = endpointExtendedOffBBX.xy + dirv1;
				gl_Position = pt;
				setBG(); EmitVertex();
				pt.xy = endpointExtendedOffBBX.xy + dirv2;
				gl_Position = pt;
				aaCutoff = 0.;
				setBG(); EmitVertex();
			} else {
				vec4 pt = endpointExtendedOffBBX;
				pt.xy = endpointExtendedOffBBX.xy - dirv2;
				gl_Position = pt;
				aaCutoff = 0.;
				setBG(); EmitVertex();
				pt.xy = endpointExtendedOffBBX.xy - dirv1;
				gl_Position = pt;
				aaCutoff = aa;
				setBG(); EmitVertex();
			}
			EndPrimitive();
			lineEdge = 0.; aaCutoff = 0.;
		} else if (connector_mode_1 > 0. && drawBackground > 0.){
			// create 1 lines, between transformedTarget -> endpointOnBBX
			// replace endpointOnBBX with midpoint of midpoints 
			vec4 midULUR = (upperLeft + upperRight) / 2.0;
			vec4 midLLLR = (lowerLeft + lowerRight) / 2.0;
			vec4 midMid = (midULUR + midLLLR) / 2.0;
			vec2 dirv;
                        drawLineAsGeometryClipped(transformedTarget, midMid, labelTop, labelBottom, labelLeft, labelRight, dirv, adjLineWidth, aa);
		} else { //connector_mode 0 or 3, or connector_mode 1 without the background
			vec2 dirv;
			drawLineAsGeometry(transformedTarget, endpointOnBBX, dirv, adjLineWidth, aa);
		}
	}
}
