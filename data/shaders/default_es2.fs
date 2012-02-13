
varying vec3 packed_data_0 ;
varying vec4 packed_data_1 ;
varying vec4 packed_data_2 ;
varying vec4 packed_data_3 ;
varying vec4 packed_data_4 ;

//varying vec3 N;
#define NORMAL packed_data_0.xyz
#define COLOR packed_data_3
#define fog (packed_data_1.w )
#define fog_color packed_data_1.x, packed_data_1.y, packed_data_1.z

uniform float fog_enabled;

uniform bool lighting_enabled;
uniform bool two_sided_lighting_enabled;
uniform bool bg_gradient;
uniform int light_count;
uniform vec4 interior_color;
uniform float interior_color_threshold;
uniform float shininess;
uniform float shininess_0;
uniform bool use_interior_color_threshold;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;
uniform int stereo_flag;
uniform mat3 matR;
uniform mat3 matL;
uniform float gamma;

vec4 ComputeColorForLight(bool is_interior, bool two_sided_lighting_enabled,
                          vec3 L, vec3 H, vec4 ambient, vec4 diffuse, float spec, float shine){
  float NdotL, NdotH;
  vec4 ret_val = vec4(0.);
  if (!is_interior){
    ret_val += ambient * COLOR;
    NdotL = dot(NORMAL, L);
    if (NdotL > 0.0) {
       ret_val += diffuse * NdotL * COLOR;
       NdotH = max(dot(NORMAL, H), 0.0);
       ret_val += spec * pow(NdotH, shine);
    }
  }
  if (two_sided_lighting_enabled && is_interior){
    NdotL = dot(-NORMAL, L);
    if (NdotL > 0.0) {
       ret_val += diffuse * NdotL * COLOR;
       NdotH = max(dot(-NORMAL, H), 0.0);
       ret_val += spec * pow(NdotH, shine);
    }
  }
  return ret_val;
}
void main()
{
  vec4 final_color = vec4(0.);

  if (lighting_enabled){
    bool is_interior = false;
    if (use_interior_color_threshold){
      vec3 viewV = vec3(0.,0.,-1.);
      float dotp = dot(NORMAL, viewV);
      is_interior = ( dotp > interior_color_threshold );
    }
    if (!two_sided_lighting_enabled && is_interior){
      final_color = interior_color;
    } else {
      final_color = (gl_LightModel.ambient) * COLOR;

      if (light_count>0){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[0].position)),
                                            normalize(vec3(gl_LightSource[0].halfVector.xyz)),
                                            gl_LightSource[0].ambient,
                                            gl_LightSource[0].diffuse,
                                            spec_value_0, shininess_0);
      if (light_count>1){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[1].position)),
                                            normalize(vec3(gl_LightSource[1].halfVector.xyz)),
                                            gl_LightSource[1].ambient,
                                            gl_LightSource[1].diffuse,
                                            spec_value, shininess);
      if (light_count>2){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[2].position)),
                                            normalize(vec3(gl_LightSource[2].halfVector.xyz)),
                                            gl_LightSource[2].ambient,
                                            gl_LightSource[2].diffuse,
                                            spec_value, shininess);
      if (light_count>3){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[3].position)),
                                            normalize(vec3(gl_LightSource[3].halfVector.xyz)),
                                            gl_LightSource[3].ambient,
                                            gl_LightSource[3].diffuse,
                                            spec_value, shininess);
      if (light_count>4){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[4].position)),
                                            normalize(vec3(gl_LightSource[4].halfVector.xyz)),
                                            gl_LightSource[4].ambient,
                                            gl_LightSource[4].diffuse,
                                            spec_value, shininess);
      if (light_count>5){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[5].position)),
                                            normalize(vec3(gl_LightSource[5].halfVector.xyz)),
                                            gl_LightSource[5].ambient,
                                            gl_LightSource[5].diffuse,
                                            spec_value, shininess);
      if (light_count>6){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[6].position)),
                                            normalize(vec3(gl_LightSource[6].halfVector.xyz)),
                                            gl_LightSource[6].ambient,
                                            gl_LightSource[6].diffuse,
                                            spec_value, shininess);
      if (light_count>7){
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled,
                                            normalize(vec3(gl_LightSource[7].position)),
                                            normalize(vec3(gl_LightSource[7].halfVector.xyz)),
                                            gl_LightSource[7].ambient,
                                            gl_LightSource[7].diffuse,
                                            spec_value, shininess);
      }}}}}}}}
    }
/*
      for (i=0; i<light_count;i++){
        vec3 L = normalize(vec3(gl_LightSource[i].position));
        vec3 H = normalize(vec3(gl_LightSource[i].halfVector.xyz));
	vec4 ambient = gl_LightSource[i].ambient;
	vec4 diffuse = gl_LightSource[i].diffuse;
        float spec = 0., shine = 0.;
        if (i==0){
          spec = spec_value_0;
	  shine = shininess_0;
        } else if (spec_count >= i){
          spec = spec_value;
	  shine = shininess;
        }
	final_color += ComputeColorForLight(is_interior, two_sided_lighting_enabled, L, H, ambient, diffuse, spec, shine);
      }

    }
*/
  } else {
    final_color = COLOR;
  }

  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  vec4 f = vec4(mix(vec3(fog_color), final_color.rgb, cfog), COLOR.a);

  if(stereo_flag==-1)
    gl_FragColor = vec4(matL * pow(f.rgb,vec3(gamma,gamma,gamma)), f.a);
  else if (stereo_flag==0)
    gl_FragColor = f;
  else if (stereo_flag==1)
    gl_FragColor = vec4(matR * pow(f.rgb, vec3(gamma,gamma,gamma)), f.a);
}

