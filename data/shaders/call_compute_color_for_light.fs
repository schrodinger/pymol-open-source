/*
 * This shader snipplet gets included multiple times (once for every light)
 *
 * postfix = "_0"       if light == 0
 * postfix = " * 0.0"   if light > spec_count
 * postfix = ""         else
 */

lighting += ComputeLighting(normal,
        g_LightSource[`light`].position.xyz,
        g_LightSource[`light`].diffuse.r,
        spec_value`postfix`,
        shininess`postfix`);
