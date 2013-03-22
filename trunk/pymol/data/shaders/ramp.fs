
varying vec4 COLOR ;
varying vec3 NORMAL ;

#include ANAGLYPH_HEADER

void main()
{
  float NdotV = dot(NORMAL, vec3(0.,0.,1.));
  vec4 fColor = vec4(NdotV * COLOR.xyz, COLOR.a);
//  vec4 fColor = NdotV * COLOR; // slightly transparent based on normal
//  vec4 fColor = vec4(NdotV * COLOR.rgb, (.5 + (.5 * NdotV)) * COLOR.a); // a little less transparent

#include ANAGLYPH_BODY
}

