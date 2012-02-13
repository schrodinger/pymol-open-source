!!ARBvp1.0

# input contains the sphere radius in model coordinates
PARAM sphereRadius = program.env[0];
PARAM half = {0.5, 0.5, 0.0, 2.0 };
PARAM zero = {0.0, 0.0, 0.0, 1.0 };

ATTRIB vertexPosition  = vertex.position;
ATTRIB vertexNormal    = vertex.normal;
ATTRIB textureCoord    = vertex.texcoord;
OUTPUT outputPosition  = result.position;

TEMP   pos, rad, shf, txt, tip;

# Transform the vertex by the modelview matrix to get into the frame of the camera

DP4    pos.x, state.matrix.modelview.row[0], vertexPosition;
DP4    pos.y, state.matrix.modelview.row[1], vertexPosition;
DP4    pos.z, state.matrix.modelview.row[2], vertexPosition;
DP4    pos.w, state.matrix.modelview.row[3], vertexPosition;

# copy current texture coords
MOV    txt.xyzw, textureCoord.xyzw;

# scale the radius by a factor of two
MUL    rad.xy, 2.0, sphereRadius.z;

# shift the texture coordinates to the origin
SUB    shf.xy, textureCoord, {0.5, 0.5, 0.0, 0.0};

# multiply them to get the vertex offset

MUL    shf.xy, rad, shf;

# define the new vertex for corner of sphere

ADD    pos.xy, pos, shf;

# apply the projection matrix to get clip coordinates 
DP4    outputPosition.x, state.matrix.projection.row[0], pos;
DP4    outputPosition.y, state.matrix.projection.row[1], pos;
DP4    shf.z, state.matrix.projection.row[2], pos;
DP4    shf.w, state.matrix.projection.row[3], pos;
MOV    outputPosition.zw, shf;

# compute camera position for front tip of the sphere
ADD    pos.z, pos.z, sphereRadius;

# compute Zc and Wc for front tip of the sphere
DP4    tip.z, state.matrix.projection.row[2], pos;
DP4    tip.w, state.matrix.projection.row[3], pos;

# compute 1/Wc for sphere tip 
RCP    rad.z, tip.w;

# put sphere center Zc into tip.w 
MOV    tip.w, shf.z;

# compute 1/Wc for sphere center 
RCP    rad.w, shf.w;

# compute Z/Wc for both sphere tip (->txt.z) and center (->txt.w) 
MUL    txt.zw, tip, rad;

# move into range 0.0-1.0 to get the normalized depth coordinate (0.5*(Zc/Wc)+0.5) 
ADD    txt.zw, {0.0,0.0,1.0,1.0}, txt;
MUL    txt.zw, {0.0,0.0,0.5,0.5}, txt;

# Pass the color through
MOV    result.color, vertex.color;

# Pass texture through
MOV    result.texcoord, txt;

END

