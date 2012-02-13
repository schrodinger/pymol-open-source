
// cylinder imposter fragment shader

uniform bool lighting_enabled;

uniform float fog_enabled;
uniform bool bg_gradient;
uniform vec3 fog_color_top;
uniform vec3 fog_color_bottom;
uniform float inv_height;
uniform float ortho;
uniform float flat_caps;
uniform bool filter_front_facing;
uniform bool two_sided_lighting_enabled;
uniform int light_count;
uniform float shininess;
uniform float shininess_0;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;
uniform float half_bond;
uniform int stereo_flag;
uniform mat3 matR;
uniform mat3 matL;
uniform float gamma;
//varying vec3 point; // surface point
//varying vec3 axis; // cylinder axis
//varying vec3 base; // cylinder base
//varying vec3 end; // cylinder end
//varying vec3 U; // cylinder base plane coordinates
//varying vec3 V;
//varying float radius; // radius
//varying float cap; // should we draw the endcap
//varying float inv_sqr_height;

varying vec4 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_3 ;
varying vec4 packed_data_4 ;
varying vec4 packed_data_5 ;

// point -> surface_point b/c preprocessor replaces _point
#define surface_point ( packed_data_0.xyz )
#define axis ( packed_data_1.xyz )
#define base ( packed_data_2.xyz )
// end -> end_cyl
#define end_cyl packed_data_3.xyz
#define U ( packed_data_4.xyz )
#define V ( packed_data_5.xyz )
#define radius ( packed_data_3.w )
#define cap ( packed_data_4.w )
#define inv_sqr_height ( packed_data_5.w )

varying vec4 color1;
varying vec4 color2;

vec4 ComputeColorForLight(vec3 N, vec3 L, vec3 H, vec4 ambient, vec4 diffuse, float spec, float shine, vec4 color){
  float NdotL, NdotH;
  vec4 ret_val = vec4(0.);
  ret_val += ambient * color;
  NdotL = dot(N, L);
  if (NdotL > 0.0) {
    ret_val += diffuse * NdotL * color;
    NdotH = max(dot(N, H), 0.0);
    ret_val += spec * pow(NdotH, shine);
  }
  return ret_val;
}

void main(void)
{
    // cull back face - otherwise we are drawing all pixels twice
    // this change gives roughly 2x speedup
    if (filter_front_facing && !gl_FrontFacing) 
      discard; 

    vec3 ray_target = surface_point;
    vec3 ray_origin = vec3(0.0);
    vec3 ray_direction = mix(normalize(ray_origin - ray_target), vec3(0.0, 0.0, 1.0), ortho);

    // basis is local system of coordinates for the cylinder
    mat3 basis = mat3(U, V, axis);

    vec3 diff = ray_target - 0.5 * (base + end_cyl);
    vec3 P = diff * basis;

    // angle (cos) between cylinder cylinder_axis and ray direction
    float dz = dot(axis, ray_direction);

    float radius2 = radius*radius;

    // calculate distance to the cylinder from ray origin
    vec3 D = vec3(dot(U, ray_direction),
                  dot(V, ray_direction),
                  dz);
    float a0 = P.x*P.x + P.y*P.y - radius2;
    float a1 = P.x*D.x + P.y*D.y;
    float a2 = D.x*D.x + D.y*D.y;
    // calculate a dicriminant of the above quadratic equation
    float d = a1*a1 - a0*a2;
    if (d < 0.0)
        // outside of the cylinder
        discard;

    float dist = (-a1 + sqrt(d))/a2;

    // point of intersection on cylinder surface
    vec3 new_point = ray_target + dist * ray_direction;

    vec3 tmp_point = new_point - base;
    vec3 normal = normalize(tmp_point - axis * dot(tmp_point, axis));

    float ratio = dot(new_point-base, vec3(end_cyl-base)) * inv_sqr_height;

    ray_origin = mix(ray_origin, surface_point, ortho);

    /* cap :  4 bits : 1st - frontcap
                       2nd - endcap
                       3rd - frontcap round
                       4th - endcap round
     */
    int icap = int(cap);
    vec4 color;
    if (icap > 15){
       float dp = clamp(-half_bond*new_point.z*inv_height, 0., .5);
       color = mix(color1, color2, smoothstep(.5 - dp, .5 + dp, ratio));
       icap = icap - 16;
    } else {
       color = mix(color1, color2, ratio);
    }
    float frontcap = 0.0, frontcapround = 0.0;
    float endcap = 0.0, endcapround = 0.0;

    if (icap == 1 || icap == 3 || icap == 5 || icap == 7 ||
        icap == 9 || icap == 11 || icap == 13 || icap == 15)
        frontcap = 1.0;

    if (icap == 2 || icap == 3 || icap == 6 || icap == 7 ||
        icap == 10 || icap == 11 || icap == 14 || icap == 15)
        endcap = 1.0;

    if (frontcap > 0.5 && flat_caps < 0.5 &&
        ((icap > 3 && icap < 8) ||
         icap > 11))
        frontcapround = 1.0;

    if (endcap > 0.5 && flat_caps < 0.5 && icap > 7)
        endcapround = 1.0;

    // test front cap
    float cap_test = dot((new_point - base), axis);

    // to calculate caps, simply check the angle between
    // the point of intersection - cylinder end vector
    // and a cap plane normal (which is the cylinder cylinder_axis)
    // if the angle < 0, the point is outside of cylinder
    // test front cap

    // flat
    if (frontcapround < 0.5 && cap_test < 0.0) {
      // ray-plane intersection
      color = color1;
      float dNV = dot(-axis, ray_direction);
      if (dNV < 0.0) discard;
      float near = dot(-axis, (base)) / dNV;
      new_point = ray_direction * near + ray_origin;
      // within the cap radius?
      if (dot(new_point - base, new_point-base) > radius2) discard;
      normal = -axis;
    }

    // round
    if (frontcapround > 0.5 && cap_test < 0.0) {
      if ( frontcap < 1.0)
        discard;
      color = color1;
      vec3 sphere_direction = mix(base, ray_origin - base, ortho);
      float b = dot(sphere_direction, ray_direction);
      float pos = b*b + radius2 -dot(sphere_direction, sphere_direction);
      if ( pos < 0.0)
        discard;
      float near = mix(b + sqrt(pos), sqrt(pos) - b, ortho);
      new_point = near * ray_direction + ray_origin;     
      normal = normalize( new_point - base ); 
    }
    // test end cap

    cap_test = dot((new_point - end_cyl), axis);

    // flat
    if (endcapround < 0.5 && cap_test > 0.0) {
      // ray-plane intersection
      color = color1;
      float dNV = dot(axis, ray_direction);
      if (dNV < 0.0) discard;
      float near = dot(axis, end_cyl) / dNV;
      new_point = ray_direction * near + ray_origin;
      // within the cap radius?
      if (dot(new_point - end_cyl, new_point-base) > radius2) discard;
      normal = axis;
    }

    // round

    if (endcapround > 0.5 && cap_test > 0.0) {
      if ( endcap < 1.0)
        discard;
      color = color2;
      vec3 sphere_direction = mix(end_cyl, ray_origin - end_cyl, ortho);
      float b = dot(sphere_direction, ray_direction);
      float pos = b*b + radius2 -dot(sphere_direction, sphere_direction);
      if ( pos < 0.0)
        discard;
      float near = mix(b + sqrt(pos), sqrt(pos) - b, ortho);
      new_point = near * ray_direction + ray_origin;
      normal = normalize( new_point - end_cyl );
    }

    vec2 clipZW = new_point.z * gl_ProjectionMatrix[2].zw +
        gl_ProjectionMatrix[3].zw;

    float depth = 0.5 + 0.5 * clipZW.x / clipZW.y;

    // this is a workaround necessary for Mac
    // otherwise the modified fragment won't clip properly

    if (depth <= 0.0)
      discard;

    if (depth >= 1.0)
      discard;

    gl_FragDepth = depth;

  vec4 final_color;
  int i;
  float NdotL, NdotH;

  final_color = (gl_LightModel.ambient) * color;

  if (light_count>0){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[0].position)),
                                        normalize(vec3(gl_LightSource[0].halfVector.xyz)),
                                        gl_LightSource[0].ambient,
                                        gl_LightSource[0].diffuse,
                                        spec_value_0, shininess_0, color);
  if (light_count>1){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[1].position)),
                                        normalize(vec3(gl_LightSource[1].halfVector.xyz)),
                                        gl_LightSource[1].ambient,
                                        gl_LightSource[1].diffuse,
                                        spec_value, shininess, color);
  if (light_count>2){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[2].position)),
                                        normalize(vec3(gl_LightSource[2].halfVector.xyz)),
                                        gl_LightSource[2].ambient,
                                        gl_LightSource[2].diffuse,
                                        spec_value, shininess, color);
  if (light_count>3){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[3].position)),
                                        normalize(vec3(gl_LightSource[3].halfVector.xyz)),
                                        gl_LightSource[3].ambient,
                                        gl_LightSource[3].diffuse,
                                        spec_value, shininess, color);
  if (light_count>4){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[4].position)),
                                        normalize(vec3(gl_LightSource[4].halfVector.xyz)),
                                        gl_LightSource[4].ambient,
                                        gl_LightSource[4].diffuse,
                                        spec_value, shininess, color);
  if (light_count>5){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[5].position)),
                                        normalize(vec3(gl_LightSource[5].halfVector.xyz)),
                                        gl_LightSource[5].ambient,
                                        gl_LightSource[5].diffuse,
                                        spec_value, shininess, color);
  if (light_count>6){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[6].position)),
                                        normalize(vec3(gl_LightSource[6].halfVector.xyz)),
                                        gl_LightSource[6].ambient,
                                        gl_LightSource[6].diffuse,
                                        spec_value, shininess, color);
  if (light_count>7){
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[7].position)),
                                        normalize(vec3(gl_LightSource[7].halfVector.xyz)),
                                        gl_LightSource[7].ambient,
                                        gl_LightSource[7].diffuse,
                                        spec_value, shininess, color);
  }}}}}}}}

/*
  for (i=0; i<light_count; i++){
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

    final_color += gl_LightSource[i].ambient * color;
    NdotL = dot(normal, L);
    if (NdotL > 0.0) {
      final_color += gl_LightSource[i].diffuse * NdotL * color;
      NdotH = max(dot(normal, H), 0.0);
      final_color += spec * pow(NdotH, shine);
    }

// I don't think we need two sided lighting for cylinders since we do not
   draw the insides of the cylinders 
    if (two_sided_lighting_enabled){
      NdotL = dot(-normal, L);
      if (NdotL > 0.0) {
        final_color += gl_LightSource[i].diffuse * NdotL * color;
        NdotH = max(dot(-normal, H), 0.0);
        final_color += spec * pow(NdotH, shine);
      }
    }
//
  }
*/
  float fog = clamp((gl_Fog.end + new_point.z) * gl_Fog.scale, 0.0, 1.0);

  fog = mix(1.0, fog, fog_enabled);
    
  vec3 fog_color;

  if (bg_gradient){
    fog_color = mix(fog_color_bottom, fog_color_top, gl_FragCoord.y * inv_height);
  } else {
    fog_color = fog_color_top;
  }

  final_color.rgb = mix(fog_color, final_color.rgb, fog);

  vec4 f = vec4(final_color.rgb, color.a);

  if(stereo_flag==-1)
    gl_FragColor = vec4(matL * pow(f.rgb,vec3(gamma,gamma,gamma)), f.a);
  else if (stereo_flag==0)
    gl_FragColor = f;
  else if (stereo_flag==1)
    gl_FragColor = vec4(matR * pow(f.rgb, vec3(gamma,gamma,gamma)), f.a);
}

