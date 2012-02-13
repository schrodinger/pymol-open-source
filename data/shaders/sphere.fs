// Sphere impostor fragment shader

uniform bool lighting_enabled;

uniform float ortho;

uniform float fog_enabled;
uniform bool bg_gradient;
uniform vec3 fog_color_top;
uniform vec3 fog_color_bottom;
uniform float inv_height;
uniform int light_count;
uniform float shininess;
uniform float shininess_0;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;
uniform int stereo_flag;
uniform mat3 matL;
uniform mat3 matR;
uniform float gamma;

varying vec4 COLOR;
varying vec3 sphere_center;
varying float radius2;
varying vec3 point;
// varying fog;

vec4 ComputeColorForLight(vec3 N, vec3 L, vec3 H, vec4 ambient, vec4 diffuse, float spec, float shine){
  float NdotL, NdotH;
  vec4 ret_val = vec4(0.);
  ret_val += ambient * COLOR;
  NdotL = dot(N, L);
  if (NdotL > 0.0) {
    ret_val += diffuse * NdotL * COLOR;
    NdotH = max(dot(N, H), 0.0);
    ret_val += spec * pow(NdotH, shine);
  }
  return ret_val;
}

void main(void)
{
    vec3 ray_origin = mix(vec3(0.,0.,0.), point, ortho);

    vec3 ray_direction = mix(normalize(point), vec3(0., 0., 1.), ortho);

    vec3 sphere_direction = mix(sphere_center, ray_origin - sphere_center, ortho);

    // Calculate sphere-ray intersection
    float b = dot(sphere_direction, ray_direction);

    float position = b * b + radius2 - dot(sphere_direction, sphere_direction);

    // Check if the ray missed the sphere
    if (position < 0.0)
       discard;

    // Calculate nearest point of intersection
    float nearest = mix(b - sqrt(position), sqrt(position) - b, ortho);

    // Calculate intersection point on the sphere surface.  The ray
    // origin is at the quad (center point), so we need to project
    // back towards the user to get the front face.
    vec3 ipoint = nearest * ray_direction + ray_origin;

    // Calculate normal at the intersection point
    vec3 N = normalize(ipoint - sphere_center);

    // Calculate depth in clipping space 
    vec2 clipZW = ipoint.z * gl_ProjectionMatrix[2].zw +
        gl_ProjectionMatrix[3].zw;

    float depth = 0.5 + 0.5 * clipZW.x / clipZW.y;

    // this is a workaround necessary for Mac
    // otherwise the modified fragment won't clip properly

    if (depth <= 0.0)
      discard;

    if (depth >= 1.0)
      discard;

    if (COLOR.a <= 1.0)
      gl_FragDepth = depth;
    else
      gl_FragDepth = 1.0;


    vec4 color;

    vec3 L0 = normalize(vec3(gl_LightSource[0].position) - ipoint);

    float NdotL = max(dot(N, L0), 0.0);
    float NdotH;

    vec4 final_color = (gl_LightModel.ambient) * COLOR;

    if (light_count>0){
      final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[0].position)),
                                          normalize(vec3(gl_LightSource[0].halfVector.xyz)),
                                          gl_LightSource[0].ambient,
                                          gl_LightSource[0].diffuse,
                                          spec_value_0, shininess_0);
      if (light_count>1){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[1].position)),
                                            normalize(vec3(gl_LightSource[1].halfVector.xyz)),
                                            gl_LightSource[1].ambient,
                                            gl_LightSource[1].diffuse,
                                            spec_value, shininess);
      if (light_count>2){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[2].position)),
                                            normalize(vec3(gl_LightSource[2].halfVector.xyz)),
                                            gl_LightSource[2].ambient,
                                            gl_LightSource[2].diffuse,
                                            spec_value, shininess);
      if (light_count>3){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[3].position)),
                                            normalize(vec3(gl_LightSource[3].halfVector.xyz)),
                                            gl_LightSource[3].ambient,
                                            gl_LightSource[3].diffuse,
                                            spec_value, shininess);
      if (light_count>4){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[4].position)),
                                            normalize(vec3(gl_LightSource[4].halfVector.xyz)),
                                            gl_LightSource[4].ambient,
                                            gl_LightSource[4].diffuse,
                                            spec_value, shininess);
      if (light_count>5){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[5].position)),
                                            normalize(vec3(gl_LightSource[5].halfVector.xyz)),
                                            gl_LightSource[5].ambient,
                                            gl_LightSource[5].diffuse,
                                            spec_value, shininess);
      if (light_count>6){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[6].position)),
                                            normalize(vec3(gl_LightSource[6].halfVector.xyz)),
                                            gl_LightSource[6].ambient,
                                            gl_LightSource[6].diffuse,
                                            spec_value, shininess);
      if (light_count>7){
	final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[7].position)),
                                            normalize(vec3(gl_LightSource[7].halfVector.xyz)),
                                            gl_LightSource[7].ambient,
                                            gl_LightSource[7].diffuse,
                                            spec_value, shininess);
    }}}}}}}}

/*
    int i;
    for (i=0; i<light_count;i++){
      vec3 L = normalize(gl_LightSource[i].position.xyz);
      vec3 H = normalize(gl_LightSource[i].halfVector.xyz);
      float spec = 0., shine = 0.;
      if (i==0){
        spec = spec_value_0;
        shine = shininess_0;
      } else if (spec_count >= i){
        spec = spec_value;
        shine = shininess;
      }
      final_color += gl_LightSource[i].ambient * COLOR;
      NdotL = dot(N, L);
      if (NdotL > 0.0) {
        final_color += gl_LightSource[i].diffuse * NdotL * COLOR;
        NdotH = max(dot(N, H), 0.0);
        final_color += spec * pow(NdotH, shine);
      }
    }
*/
    float fog = clamp((gl_Fog.end + ipoint.z) * gl_Fog.scale, 0.0, 1.0);
    fog = mix(1.0, fog, fog_enabled);
    vec3 fog_color;

    if (bg_gradient){
      fog_color = mix(fog_color_bottom, fog_color_top, gl_FragCoord.y * inv_height);
    } else {
      fog_color = fog_color_top;
    }

    final_color.rgb = mix(fog_color, final_color.rgb, fog);

    vec4 f = vec4(final_color.rgb, COLOR.a);

  if(stereo_flag==-1)
    gl_FragColor = vec4(matL * pow(f.rgb,vec3(gamma,gamma,gamma)), f.a);
  else if (stereo_flag==0)
    gl_FragColor = f;
  else if (stereo_flag==1)
    gl_FragColor = vec4(matR * pow(f.rgb, vec3(gamma,gamma,gamma)), f.a);
}
