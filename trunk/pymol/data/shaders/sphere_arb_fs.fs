!!ARBfp1.0

PARAM fogInfo = program.env[0];
PARAM fogColor = state.fog.color;
ATTRIB fogCoord = fragment.fogcoord;

TEMP pln, norm, depth, color, light, spec, fogFactor;

# fully clip spheres that hit the camera
KIL fragment.texcoord.z;

# move texture coordinates to origin

MOV norm.z, 0;
SUB norm.xy, fragment.texcoord, {0.5,0.5,0.0,0.0};

# compute x^2 + y^2, if > 0.25 then kill the pixel -- not in sphere

# kill pixels that aren't in the center circle
DP3 pln.z, norm, norm;
SUB pln.z, 0.25, pln.z;
KIL pln.z;

# build a complete unit normal
MUL pln.z, 4.0, pln.z;
RSQ pln.z, pln.z;
MUL norm.xy, 2.0, norm;
RCP norm.z, pln.z;

# interpolate the Zndc coordinate on the sphere 
LRP depth.z, norm.z, fragment.texcoord.z, fragment.texcoord.w;
MOV result.depth.z, depth.z;

# light0

DP3 light, state.light[1].half, norm;
MOV light.w, 60.0;
LIT light, light;

# ambient
MOV color.xyzw, {0.06,0.06,0.06,1.0};
ADD color.xyz, light.y, 0.1;
MUL color.xyz, fragment.color, color;
MUL spec.xyz, light.z, 0.5;
ADD color.xyz, color,spec;

# apply fog using linear interp over Zndc
MAX fogFactor.x, depth.z, fogInfo.x;
SUB fogFactor.x, fogFactor.x, fogInfo.x;
MUL fogFactor.x, fogFactor.x, fogInfo.y;
LRP color.xyz, fogFactor.x, fogColor, color;
MOV result.color, color;

END

