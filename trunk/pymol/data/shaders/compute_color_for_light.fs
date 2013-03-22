#ifdef default
vec4 ComputeColorForLight(bool is_interior,
                          vec3 L, vec3 H, vec4 ambient, vec4 diffuse, float spec, float shine){
#endif
#ifdef sphere
vec4 ComputeColorForLight(vec3 NORMAL, vec3 L, vec3 H, vec4 ambient, vec4 diffuse, float spec, float shine){
#endif
#ifdef cylinder
vec4 ComputeColorForLight(vec3 NORMAL, vec3 L, vec3 H, vec4 ambient, vec4 diffuse, float spec, float shine, vec4 COLOR){
#endif
  float NdotL, NdotH;
  vec4 ret_val = vec4(0.);
#ifdef default
  if (!is_interior){
#endif
    ret_val += ambient * COLOR;
    NdotL = dot(NORMAL, L);
    if (NdotL > 0.0) {
       ret_val += diffuse * NdotL * COLOR;
       NdotH = max(dot(NORMAL, H), 0.0);
       ret_val += spec * pow(NdotH, shine);
    }
#ifdef default
  }
  if (two_sided_lighting_enabled && is_interior){
    NdotL = dot(-NORMAL, L);
    if (NdotL > 0.0) {
       ret_val += diffuse * NdotL * COLOR;
       NdotH = max(dot(-NORMAL, H), 0.0);
       ret_val += spec * pow(NdotH, shine);
    }
  }
#endif
  return ret_val;
}
