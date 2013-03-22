
// Sphere impostor vertex shader

uniform bool lighting_enabled;
uniform float sphere_size_scale;

uniform float horizontal_adjustment;
uniform float vertical_adjustment;

attribute vec4 a_vertex_radius;
attribute vec4 a_Color;
attribute float a_rightUpFlags;

varying vec4 COLOR;
varying vec3 sphere_center;
varying float radius2;
varying vec3 point;
varying vec2 bgTextureLookup;

void main(void)
{
    // Get billboard attributes
    float radius = a_vertex_radius.w * sphere_size_scale;
    float right = -1. + 2.*mod(a_rightUpFlags, 2.);
    float up = -1. + 2.*floor(mod(a_rightUpFlags/2., 2.));
    vec4 a_Vertex = vec4(a_vertex_radius.xyz, 1.);

    COLOR = a_Color;
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
    vec3 corner_direction =  (vertical_adjustment*up) * up_vector + (horizontal_adjustment*right) * right_vector;

    // Calculate vertex of screen-oriented quad (billboard)
    vec4 vertex = vec4(a_Vertex.xyz + radius * corner_direction, 1.);

    // Calculate vertex position in modelview space
    vec4 eye_space_pos = gl_ModelViewMatrix * vertex;

    // Compute sphere position in modelview space
    vec4 tmppos = gl_ModelViewMatrix * a_Vertex;
    sphere_center = tmppos.xyz / tmppos.w;

    // Compute ray direction and origin point
    point = eye_space_pos.xyz / eye_space_pos.w;

    // Pass the transformed vertex for clipping plane calculations
    gl_ClipVertex = eye_space_pos;

    // Pass transformed vertex
    gl_Position = gl_ModelViewProjectionMatrix * vertex;
    bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}

