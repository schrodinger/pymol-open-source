// if light is 0, postfix = _0, otherwise blank
#ifdef sphere
    final_color += ComputeColorForLight(N, normalize(vec3(gl_LightSource[`light`].position)),
                                        normalize(vec3(gl_LightSource[`light`].halfVector.xyz)),
                                        gl_LightSource[`light`].ambient,
                                        gl_LightSource[`light`].diffuse,
                                        spec_value`postfix`, shininess`postfix`);
#endif
#ifdef default
    final_color += ComputeColorForLight(is_interior,
                                        normalize(vec3(gl_LightSource[`light`].position)),
                                        normalize(vec3(gl_LightSource[`light`].halfVector.xyz)),
                                        gl_LightSource[`light`].ambient,
                                        gl_LightSource[`light`].diffuse,
                                        spec_value`postfix`, shininess`postfix`);
#endif
#ifdef cylinder
    final_color += ComputeColorForLight(normal, normalize(vec3(gl_LightSource[`light`].position)),
                                        normalize(vec3(gl_LightSource[`light`].halfVector.xyz)),
                                        gl_LightSource[`light`].ambient,
                                        gl_LightSource[`light`].diffuse,
                                        spec_value`postfix`, shininess`postfix`, color);
#endif
