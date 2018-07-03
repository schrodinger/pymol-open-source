#include webgl_header.vs

// cylinder imposter vertex shader

attribute vec3 attr_vertex1;
attribute vec3 attr_vertex2;
attribute vec4 a_Color;
attribute vec4 a_Color2;
attribute float attr_radius;
attribute float a_cap;

attribute float attr_flags;
varying vec3 surface_point ;
varying vec3 axis ;
varying vec3 base ;
varying vec3 end_cyl ;
varying vec3 U ;
varying vec3 V ;
varying float radius;
varying float cap;
varying float inv_sqr_height;
varying vec4 color1;
varying vec4 color2;
varying vec2 bgTextureLookup;

uniform float uni_radius;

// get_bit_and_shift: returns 0 or 1
float get_bit_and_shift(inout float bits) {
  float bit = mod(bits, 2.0);
  bits = (bits - bit) / 2.0;
  return step(.5, bit);
}

void main(void)
{
    float uniformglscale = length(g_NormalMatrix[0]);

    if (uni_radius!=0.0){
        radius = uni_radius * attr_radius;
    } else {
        radius = attr_radius;
    }

    color1 = a_Color;
    color2 = a_Color2;

    vec3 attr_axis = attr_vertex2 - attr_vertex1;

    cap = a_cap;

    // calculate reciprocal of squared height
    inv_sqr_height = length(attr_axis) / uniformglscale;
    inv_sqr_height *= inv_sqr_height;
    inv_sqr_height = 1.0 / inv_sqr_height;

    // h is a normalized cylinder axis
    vec3 h = normalize(attr_axis);
    // axis is the cylinder axis in modelview coordinates
    axis = normalize(g_NormalMatrix * h);
    // u, v, h is local system of coordinates
    vec3 u = cross(h, vec3(1.0, 0.0, 0.0));
    if (dot(u,u) < 0.001) 
      u = cross(h, vec3(0.0, 1.0, 0.0));
    u = normalize(u);
    vec3 v = normalize(cross(u, h));

    // transform to modelview coordinates
    U = normalize(g_NormalMatrix * u);
    V = normalize(g_NormalMatrix * v);

    vec4 base4 = g_ModelViewMatrix * vec4(attr_vertex1, 1.0);
    base = base4.xyz;
    vec4 end4 = g_ModelViewMatrix * vec4(attr_vertex2, 1.0);
    end_cyl = end4.xyz;

    // compute bounding box vertex position
    vec4 vertex = vec4(attr_vertex1, 1.0); 
    float packed_flags = attr_flags;
    float out_v = get_bit_and_shift(packed_flags);
    float up_v = get_bit_and_shift(packed_flags);
    float right_v = get_bit_and_shift(packed_flags);
    vertex.xyz += up_v * attr_axis;
    vertex.xyz += (2.0 * right_v - 1.0) * radius * u;
    vertex.xyz += (2.0 * out_v - 1.0) * radius * v;
    vertex.xyz += (2.0 * up_v - 1.0) * radius * h;

    vec4 tvertex = g_ModelViewMatrix * vertex;
    surface_point = tvertex.xyz;

    gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * vertex;

    // support uniform scaling
    radius /= uniformglscale;

    // clamp z on front clipping plane if impostor box would be clipped.
    // (we ultimatly want to clip on the calculated depth in the fragment
    // shader, not the depth of the box face)
    if (gl_Position.z / gl_Position.w < -1.0) {
        // upper bound of possible cylinder z extend
        float diff = abs(base4.z - end4.z) + radius * 3.5;

        // z-`diff`-offsetted vertex
        vec4 inset = g_ModelViewMatrix * vertex;
        inset.z -= diff;
        inset = g_ProjectionMatrix * inset;

        // if offsetted vertex is within front clipping plane, then clamp
        if (inset.z / inset.w > -1.0) {
            gl_Position.z = -gl_Position.w;
        }
    }
    bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;

}

