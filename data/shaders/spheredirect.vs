
// Sphere impostor vertex shader

uniform bool lighting_enabled;
uniform float sphere_size_scale;

attribute vec4 sphere_attributes;

varying vec4 COLOR;
varying vec3 sphere_center;
varying float radius2;
varying vec3 point;

void main(void)
{
    // Get billboard attributes

    float right = sphere_attributes.x;
    float up = sphere_attributes.y;
    float radius = sphere_attributes.z * sphere_size_scale;
    COLOR = gl_Color;
    radius2 = radius * radius; // compute squared radius 

    vec3 right_vector = vec3(gl_ModelViewMatrix[0][0],
            gl_ModelViewMatrix[1][0], gl_ModelViewMatrix[2][0]);

    vec3 up_vector = vec3(gl_ModelViewMatrix[0][1],
            gl_ModelViewMatrix[1][1], gl_ModelViewMatrix[2][1]);

    // We need to project the vertex out to the edge of the square, which
    // is the following distance:
    // float corner_distance = sqrt(2.0 * radius2);
    // but since we need to normalize the corner vector computed below
    // which has length sqrt(2.0), we can simply use radius as corner distance
    // to compute vertex position of screen-oriented quad.

    // Compute corner vector
    vec3 corner_direction =  up * up_vector + right * right_vector;

    // Calculate vertex of screen-oriented quad (billboard)
    vec4 vertex = vec4(gl_Vertex.xyz + radius * corner_direction, 1.);

    // Calculate vertex position in modelview space
    vec4 eye_space_pos = gl_ModelViewMatrix * vertex;

    // Compute sphere position in modelview space
    vec4 tmppos = gl_ModelViewMatrix * gl_Vertex;
    sphere_center = vec3(tmppos) / tmppos.w;

    // Compute ray direction and origin point
    point = vec3(eye_space_pos) / eye_space_pos.w;

    // Pass fog coordinate
    gl_FogFragCoord = abs(sphere_center.z); 

    // Pass the transformed vertex for clipping plane calculations
    gl_ClipVertex = eye_space_pos;

    vec3 eye_pos = vec3(gl_ModelViewMatrix * gl_Vertex);

    // Pass transformed vertex
    gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * vertex;
}

