
/* BEGIN GLUT EXCERPT.  THE FOLLOWING CODE IS:

 * Copyright (c) Mark J. Kilgard, 1994. 
 
 * This program is freely distributable without licensing fees 
 * and is provided without guarantee or warrantee expressed or 
 * implied. This program is -not- in the public domain. 

 * Modifications by Warren L. DeLano, 2004.  

*/

#include "FontGLUT.h"

/* GENERATED FILE -- DO NOT MODIFY */

#ifdef _WIN32
/* XXX Work around Microsoft OpenGL 1.1 bug where glBitmap with
   a height or width of zero does not advance the raster position
   as specified by OpenGL. (Cosmo OpenGL does not have this bug.) */
static const unsigned char  ch0data[] = { 0x0 };
static const FontGLUTBitmapCharRec ch0 = {1,1,0,0,8,ch0data};
#else
static const FontGLUTBitmapCharRec ch0 = {0,0,0,0,8,0};
#endif

#ifdef _WIN32
/* XXX Work around Microsoft OpenGL 1.1 bug where glBitmap with
   a height or width of zero does not advance the raster position
   as specified by OpenGL. (Cosmo OpenGL does not have this bug.) */
static const unsigned char  ch32data[] = { 0x0 };
static const FontGLUTBitmapCharRec ch32 = {1,1,0,0,8,ch32data};
#else
static const FontGLUTBitmapCharRec ch32 = {0,0,0,0,8,0};
#endif

#ifdef _WIN32
/* XXX Work around Microsoft OpenGL 1.1 bug where glBitmap with
   a height or width of zero does not advance the raster position
   as specified by OpenGL. (Cosmo OpenGL does not have this bug.) */
static const unsigned char  ch127data[] = { 0x0 };
static const FontGLUTBitmapCharRec ch127 = {1,1,0,0,8,ch127data};
#else
static const FontGLUTBitmapCharRec ch127 = {0,0,0,0,8,0};
#endif

#ifdef _WIN32
/* XXX Work around Microsoft OpenGL 1.1 bug where glBitmap with
   a height or width of zero does not advance the raster position
   as specified by OpenGL. (Cosmo OpenGL does not have this bug.) */
static const unsigned char  ch160data[] = { 0x0 };
static const FontGLUTBitmapCharRec ch160 = {1,1,0,0,8,ch160data};
#else
static const FontGLUTBitmapCharRec ch160 = {0,0,0,0,8,0};
#endif

/* char: 0xff */

static const unsigned char  ch255data[] = {
0x78,0x84,0x4,0x74,0x8c,0x84,0x84,0x84,0x0,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch255 = {6,12,-1,2,8,ch255data};

/* char: 0xfe */

static const unsigned char  ch254data[] = {
0x80,0x80,0xb8,0xc4,0x84,0x84,0xc4,0xb8,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch254 = {6,10,-1,2,8,ch254data};

/* char: 0xfd */

static const unsigned char  ch253data[] = {
0x78,0x84,0x4,0x74,0x8c,0x84,0x84,0x84,0x0,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch253 = {6,12,-1,2,8,ch253data};

/* char: 0xfc */

static const unsigned char  ch252data[] = {
0x74,0x88,0x88,0x88,0x88,0x88,0x0,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch252 = {6,10,-1,0,8,ch252data};

/* char: 0xfb */

static const unsigned char  ch251data[] = {
0x74,0x88,0x88,0x88,0x88,0x88,0x0,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch251 = {6,10,-1,0,8,ch251data};

/* char: 0xfa */

static const unsigned char  ch250data[] = {
0x74,0x88,0x88,0x88,0x88,0x88,0x0,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch250 = {6,10,-1,0,8,ch250data};

/* char: 0xf9 */

static const unsigned char  ch249data[] = {
0x74,0x88,0x88,0x88,0x88,0x88,0x0,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch249 = {6,10,-1,0,8,ch249data};

/* char: 0xf8 */

static const unsigned char  ch248data[] = {
0x80,0x78,0xc4,0xa4,0x94,0x8c,0x78,0x4,
};

static const FontGLUTBitmapCharRec ch248 = {6,8,-1,1,8,ch248data};

/* char: 0xf7 */

static const unsigned char  ch247data[] = {
0x20,0x20,0x0,0xf8,0x0,0x20,0x20,
};

static const FontGLUTBitmapCharRec ch247 = {5,7,-1,-1,8,ch247data};

/* char: 0xf6 */

static const unsigned char  ch246data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,0x0,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch246 = {6,10,-1,0,8,ch246data};

/* char: 0xf5 */

static const unsigned char  ch245data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,0x0,0x0,0x50,0x28,
};

static const FontGLUTBitmapCharRec ch245 = {6,10,-1,0,8,ch245data};

/* char: 0xf4 */

static const unsigned char  ch244data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,0x0,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch244 = {6,10,-1,0,8,ch244data};

/* char: 0xf3 */

static const unsigned char  ch243data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,0x0,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch243 = {6,10,-1,0,8,ch243data};

/* char: 0xf2 */

static const unsigned char  ch242data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,0x0,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch242 = {6,10,-1,0,8,ch242data};

/* char: 0xf1 */

static const unsigned char  ch241data[] = {
0x84,0x84,0x84,0x84,0xc4,0xb8,0x0,0x0,0x50,0x28,
};

static const FontGLUTBitmapCharRec ch241 = {6,10,-1,0,8,ch241data};

/* char: 0xf0 */

static const unsigned char  ch240data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,0x8,0x50,0x30,0x48,
};

static const FontGLUTBitmapCharRec ch240 = {6,10,-1,0,8,ch240data};

/* char: 0xef */

static const unsigned char  ch239data[] = {
0xf8,0x20,0x20,0x20,0x20,0x60,0x0,0x0,0x50,0x50,
};

static const FontGLUTBitmapCharRec ch239 = {5,10,-1,0,8,ch239data};

/* char: 0xee */

static const unsigned char  ch238data[] = {
0xf8,0x20,0x20,0x20,0x20,0x60,0x0,0x0,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch238 = {5,10,-1,0,8,ch238data};

/* char: 0xed */

static const unsigned char  ch237data[] = {
0xf8,0x20,0x20,0x20,0x20,0x60,0x0,0x0,0x40,0x20,
};

static const FontGLUTBitmapCharRec ch237 = {5,10,-1,0,8,ch237data};

/* char: 0xec */

static const unsigned char  ch236data[] = {
0xf8,0x20,0x20,0x20,0x20,0x60,0x0,0x0,0x20,0x40,
};

static const FontGLUTBitmapCharRec ch236 = {5,10,-1,0,8,ch236data};

/* char: 0xeb */

static const unsigned char  ch235data[] = {
0x78,0x84,0x80,0xfc,0x84,0x78,0x0,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch235 = {6,10,-1,0,8,ch235data};

/* char: 0xea */

static const unsigned char  ch234data[] = {
0x78,0x84,0x80,0xfc,0x84,0x78,0x0,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch234 = {6,10,-1,0,8,ch234data};

/* char: 0xe9 */

static const unsigned char  ch233data[] = {
0x78,0x84,0x80,0xfc,0x84,0x78,0x0,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch233 = {6,10,-1,0,8,ch233data};

/* char: 0xe8 */

static const unsigned char  ch232data[] = {
0x78,0x84,0x80,0xfc,0x84,0x78,0x0,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch232 = {6,10,-1,0,8,ch232data};

/* char: 0xe7 */

static const unsigned char  ch231data[] = {
0x20,0x10,0x78,0x84,0x80,0x80,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch231 = {6,8,-1,2,8,ch231data};

/* char: 0xe6 */

static const unsigned char  ch230data[] = {
0x6c,0x92,0x90,0x7c,0x12,0x6c,
};

static const FontGLUTBitmapCharRec ch230 = {7,6,0,0,8,ch230data};

/* char: 0xe5 */

static const unsigned char  ch229data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,0x0,0x30,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch229 = {6,10,-1,0,8,ch229data};

/* char: 0xe4 */

static const unsigned char  ch228data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,0x0,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch228 = {6,10,-1,0,8,ch228data};

/* char: 0xe3 */

static const unsigned char  ch227data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,0x0,0x0,0x50,0x28,
};

static const FontGLUTBitmapCharRec ch227 = {6,10,-1,0,8,ch227data};

/* char: 0xe2 */

static const unsigned char  ch226data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,0x0,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch226 = {6,10,-1,0,8,ch226data};

/* char: 0xe1 */

static const unsigned char  ch225data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,0x0,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch225 = {6,10,-1,0,8,ch225data};

/* char: 0xe0 */

static const unsigned char  ch224data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,0x0,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch224 = {6,10,-1,0,8,ch224data};

/* char: 0xdf */

static const unsigned char  ch223data[] = {
0x80,0xb8,0xc4,0x84,0x84,0xf8,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch223 = {6,9,-1,1,8,ch223data};

/* char: 0xde */

static const unsigned char  ch222data[] = {
0x80,0x80,0x80,0xf8,0x84,0x84,0x84,0xf8,0x80,
};

static const FontGLUTBitmapCharRec ch222 = {6,9,-1,0,8,ch222data};

/* char: 0xdd */

static const unsigned char  ch221data[] = {
0x20,0x20,0x20,0x20,0x50,0x88,0x88,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch221 = {5,10,-1,0,8,ch221data};

/* char: 0xdc */

static const unsigned char  ch220data[] = {
0x78,0x84,0x84,0x84,0x84,0x84,0x84,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch220 = {6,10,-1,0,8,ch220data};

/* char: 0xdb */

static const unsigned char  ch219data[] = {
0x78,0x84,0x84,0x84,0x84,0x84,0x84,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch219 = {6,10,-1,0,8,ch219data};

/* char: 0xda */

static const unsigned char  ch218data[] = {
0x78,0x84,0x84,0x84,0x84,0x84,0x84,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch218 = {6,10,-1,0,8,ch218data};

/* char: 0xd9 */

static const unsigned char  ch217data[] = {
0x78,0x84,0x84,0x84,0x84,0x84,0x84,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch217 = {6,10,-1,0,8,ch217data};

/* char: 0xd8 */

static const unsigned char  ch216data[] = {
0x80,0x78,0xc4,0xa4,0xa4,0xa4,0x94,0x94,0x8c,0x78,0x4,
};

static const FontGLUTBitmapCharRec ch216 = {6,11,-1,1,8,ch216data};

/* char: 0xd7 */

static const unsigned char  ch215data[] = {
0x84,0x48,0x30,0x30,0x48,0x84,
};

static const FontGLUTBitmapCharRec ch215 = {6,6,-1,-1,8,ch215data};

/* char: 0xd6 */

static const unsigned char  ch214data[] = {
0x7c,0x82,0x82,0x82,0x82,0x82,0x7c,0x0,0x28,0x28,
};

static const FontGLUTBitmapCharRec ch214 = {7,10,0,0,8,ch214data};

/* char: 0xd5 */

static const unsigned char  ch213data[] = {
0x7c,0x82,0x82,0x82,0x82,0x82,0x7c,0x0,0x28,0x14,
};

static const FontGLUTBitmapCharRec ch213 = {7,10,0,0,8,ch213data};

/* char: 0xd4 */

static const unsigned char  ch212data[] = {
0x7c,0x82,0x82,0x82,0x82,0x82,0x7c,0x0,0x24,0x18,
};

static const FontGLUTBitmapCharRec ch212 = {7,10,0,0,8,ch212data};

/* char: 0xd3 */

static const unsigned char  ch211data[] = {
0x7c,0x82,0x82,0x82,0x82,0x82,0x7c,0x0,0x10,0x8,
};

static const FontGLUTBitmapCharRec ch211 = {7,10,0,0,8,ch211data};

/* char: 0xd2 */

static const unsigned char  ch210data[] = {
0x7c,0x82,0x82,0x82,0x82,0x82,0x7c,0x0,0x8,0x10,
};

static const FontGLUTBitmapCharRec ch210 = {7,10,0,0,8,ch210data};

/* char: 0xd1 */

static const unsigned char  ch209data[] = {
0x82,0x86,0x8a,0x92,0xa2,0xc2,0x82,0x0,0x28,0x14,
};

static const FontGLUTBitmapCharRec ch209 = {7,10,0,0,8,ch209data};

/* char: 0xd0 */

static const unsigned char  ch208data[] = {
0xfc,0x42,0x42,0x42,0xe2,0x42,0x42,0x42,0xfc,
};

static const FontGLUTBitmapCharRec ch208 = {7,9,0,0,8,ch208data};

/* char: 0xcf */

static const unsigned char  ch207data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0xf8,0x0,0x50,0x50,
};

static const FontGLUTBitmapCharRec ch207 = {5,10,-1,0,8,ch207data};

/* char: 0xce */

static const unsigned char  ch206data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0xf8,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch206 = {5,10,-1,0,8,ch206data};

/* char: 0xcd */

static const unsigned char  ch205data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0xf8,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch205 = {5,10,-1,0,8,ch205data};

/* char: 0xcc */

static const unsigned char  ch204data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0xf8,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch204 = {5,10,-1,0,8,ch204data};

/* char: 0xcb */

static const unsigned char  ch203data[] = {
0xfc,0x80,0x80,0xf0,0x80,0x80,0xfc,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch203 = {6,10,-1,0,8,ch203data};

/* char: 0xca */

static const unsigned char  ch202data[] = {
0xfc,0x80,0x80,0xf0,0x80,0x80,0xfc,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch202 = {6,10,-1,0,8,ch202data};

/* char: 0xc9 */

static const unsigned char  ch201data[] = {
0xfc,0x80,0x80,0xf0,0x80,0x80,0xfc,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch201 = {6,10,-1,0,8,ch201data};

/* char: 0xc8 */

static const unsigned char  ch200data[] = {
0xfc,0x80,0x80,0xf0,0x80,0x80,0xfc,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch200 = {6,10,-1,0,8,ch200data};

/* char: 0xc7 */

static const unsigned char  ch199data[] = {
0x20,0x10,0x78,0x84,0x80,0x80,0x80,0x80,0x80,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch199 = {6,11,-1,2,8,ch199data};

/* char: 0xc6 */

static const unsigned char  ch198data[] = {
0x9e,0x90,0x90,0xf0,0x9c,0x90,0x90,0x90,0x6e,
};

static const FontGLUTBitmapCharRec ch198 = {7,9,0,0,8,ch198data};

/* char: 0xc5 */

static const unsigned char  ch197data[] = {
0x84,0x84,0xfc,0x84,0x84,0x48,0x30,0x30,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch197 = {6,10,-1,0,8,ch197data};

/* char: 0xc4 */

static const unsigned char  ch196data[] = {
0x84,0x84,0xfc,0x84,0x84,0x48,0x30,0x0,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch196 = {6,10,-1,0,8,ch196data};

/* char: 0xc3 */

static const unsigned char  ch195data[] = {
0x84,0x84,0xfc,0x84,0x84,0x48,0x30,0x0,0x50,0x28,
};

static const FontGLUTBitmapCharRec ch195 = {6,10,-1,0,8,ch195data};

/* char: 0xc2 */

static const unsigned char  ch194data[] = {
0x84,0x84,0xfc,0x84,0x84,0x48,0x30,0x0,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch194 = {6,10,-1,0,8,ch194data};

/* char: 0xc1 */

static const unsigned char  ch193data[] = {
0x84,0x84,0xfc,0x84,0x84,0x48,0x30,0x0,0x20,0x10,
};

static const FontGLUTBitmapCharRec ch193 = {6,10,-1,0,8,ch193data};

/* char: 0xc0 */

static const unsigned char  ch192data[] = {
0x84,0x84,0xfc,0x84,0x84,0x48,0x30,0x0,0x10,0x20,
};

static const FontGLUTBitmapCharRec ch192 = {6,10,-1,0,8,ch192data};

/* char: 0xbf */

static const unsigned char  ch191data[] = {
0x78,0x84,0x84,0x80,0x40,0x20,0x20,0x0,0x20,
};

static const FontGLUTBitmapCharRec ch191 = {6,9,-1,0,8,ch191data};

/* char: 0xbe */

static const unsigned char  ch190data[] = {
0x6,0x1a,0x12,0xa,0x66,0x92,0x10,0x20,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch190 = {7,10,0,0,8,ch190data};

/* char: 0xbd */

static const unsigned char  ch189data[] = {
0x1e,0x10,0xc,0x2,0xf2,0x4c,0x40,0x40,0xc0,0x40,
};

static const FontGLUTBitmapCharRec ch189 = {7,10,0,0,8,ch189data};

/* char: 0xbc */

static const unsigned char  ch188data[] = {
0x6,0x1a,0x12,0xa,0xe6,0x42,0x40,0x40,0xc0,0x40,
};

static const FontGLUTBitmapCharRec ch188 = {7,10,0,0,8,ch188data};

/* char: 0xbb */

static const unsigned char  ch187data[] = {
0x90,0x48,0x24,0x12,0x24,0x48,0x90,
};

static const FontGLUTBitmapCharRec ch187 = {7,7,0,-1,8,ch187data};

/* char: 0xba */

static const unsigned char  ch186data[] = {
0xf0,0x0,0x60,0x90,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch186 = {4,6,-1,-3,8,ch186data};

/* char: 0xb9 */

static const unsigned char  ch185data[] = {
0xe0,0x40,0x40,0x40,0xc0,0x40,
};

static const FontGLUTBitmapCharRec ch185 = {3,6,-1,-4,8,ch185data};

/* char: 0xb8 */

static const unsigned char  ch184data[] = {
0xc0,0x40,
};

static const FontGLUTBitmapCharRec ch184 = {2,2,-3,2,8,ch184data};

/* char: 0xb7 */

static const unsigned char  ch183data[] = {
0xc0,
};

static const FontGLUTBitmapCharRec ch183 = {2,1,-3,-4,8,ch183data};

/* char: 0xb6 */

static const unsigned char  ch182data[] = {
0x28,0x28,0x28,0x28,0x68,0xe8,0xe8,0xe8,0x7c,
};

static const FontGLUTBitmapCharRec ch182 = {6,9,-1,0,8,ch182data};

/* char: 0xb5 */

static const unsigned char  ch181data[] = {
0x80,0xb4,0xcc,0x84,0x84,0x84,0x84,
};

static const FontGLUTBitmapCharRec ch181 = {6,7,-1,1,8,ch181data};

/* char: 0xb4 */

static const unsigned char  ch180data[] = {
0x80,0x40,
};

static const FontGLUTBitmapCharRec ch180 = {2,2,-3,-8,8,ch180data};

/* char: 0xb3 */

static const unsigned char  ch179data[] = {
0x60,0x90,0x10,0x20,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch179 = {4,6,-1,-4,8,ch179data};

/* char: 0xb2 */

static const unsigned char  ch178data[] = {
0xf0,0x80,0x60,0x10,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch178 = {4,6,-1,-4,8,ch178data};

/* char: 0xb1 */

static const unsigned char  ch177data[] = {
0xf8,0x0,0x20,0x20,0xf8,0x20,0x20,
};

static const FontGLUTBitmapCharRec ch177 = {5,7,-1,-1,8,ch177data};

/* char: 0xb0 */

static const unsigned char  ch176data[] = {
0x60,0x90,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch176 = {4,4,-2,-5,8,ch176data};

/* char: 0xaf */

static const unsigned char  ch175data[] = {
0xfc,
};

static const FontGLUTBitmapCharRec ch175 = {6,1,-1,-8,8,ch175data};

/* char: 0xae */

static const unsigned char  ch174data[] = {
0x38,0x44,0xaa,0xb2,0xaa,0xaa,0x92,0x44,0x38,
};

static const FontGLUTBitmapCharRec ch174 = {7,9,0,-1,8,ch174data};

/* char: 0xad */

static const unsigned char  ch173data[] = {
0xfc,
};

static const FontGLUTBitmapCharRec ch173 = {6,1,-1,-4,8,ch173data};

/* char: 0xac */

static const unsigned char  ch172data[] = {
0x4,0x4,0x4,0xfc,
};

static const FontGLUTBitmapCharRec ch172 = {6,4,-1,-1,8,ch172data};

/* char: 0xab */

static const unsigned char  ch171data[] = {
0x12,0x24,0x48,0x90,0x48,0x24,0x12,
};

static const FontGLUTBitmapCharRec ch171 = {7,7,0,-1,8,ch171data};

/* char: 0xaa */

static const unsigned char  ch170data[] = {
0xf8,0x0,0x78,0x88,0x78,0x8,0x70,
};

static const FontGLUTBitmapCharRec ch170 = {5,7,-1,-2,8,ch170data};

/* char: 0xa9 */

static const unsigned char  ch169data[] = {
0x38,0x44,0x92,0xaa,0xa2,0xaa,0x92,0x44,0x38,
};

static const FontGLUTBitmapCharRec ch169 = {7,9,0,-1,8,ch169data};

/* char: 0xa8 */

static const unsigned char  ch168data[] = {
0xd8,
};

static const FontGLUTBitmapCharRec ch168 = {5,1,-1,-8,8,ch168data};

/* char: 0xa7 */

static const unsigned char  ch167data[] = {
0x60,0x90,0x10,0x60,0x90,0x90,0x60,0x80,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch167 = {4,10,-2,0,8,ch167data};

/* char: 0xa6 */

static const unsigned char  ch166data[] = {
0x80,0x80,0x80,0x80,0x0,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch166 = {1,9,-3,0,8,ch166data};

/* char: 0xa5 */

static const unsigned char  ch165data[] = {
0x10,0x10,0x7c,0x10,0x7c,0x28,0x44,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch165 = {7,9,0,0,8,ch165data};

/* char: 0xa4 */

static const unsigned char  ch164data[] = {
0x84,0x78,0x48,0x48,0x78,0x84,
};

static const FontGLUTBitmapCharRec ch164 = {6,6,-1,-1,8,ch164data};

/* char: 0xa3 */

static const unsigned char  ch163data[] = {
0xdc,0x62,0x20,0x20,0x20,0x70,0x20,0x22,0x1c,
};

static const FontGLUTBitmapCharRec ch163 = {7,9,0,0,8,ch163data};

/* char: 0xa2 */

static const unsigned char  ch162data[] = {
0x20,0x70,0xa8,0xa0,0xa0,0xa8,0x70,0x20,
};

static const FontGLUTBitmapCharRec ch162 = {5,8,-1,-1,8,ch162data};

/* char: 0xa1 */

static const unsigned char  ch161data[] = {
0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x0,0x80,
};

static const FontGLUTBitmapCharRec ch161 = {1,9,-3,0,8,ch161data};

/* char: 0x7e '~' */

static const unsigned char  ch126data[] = {
0x90,0xa8,0x48,
};

static const FontGLUTBitmapCharRec ch126 = {5,3,-1,-6,8,ch126data};

/* char: 0x7d '}' */

static const unsigned char  ch125data[] = {
0xe0,0x10,0x10,0x20,0x18,0x20,0x10,0x10,0xe0,
};

static const FontGLUTBitmapCharRec ch125 = {5,9,-1,0,8,ch125data};

/* char: 0x7c '|' */

static const unsigned char  ch124data[] = {
0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch124 = {1,9,-3,0,8,ch124data};

/* char: 0x7b '{' */

static const unsigned char  ch123data[] = {
0x38,0x40,0x40,0x20,0xc0,0x20,0x40,0x40,0x38,
};

static const FontGLUTBitmapCharRec ch123 = {5,9,-2,0,8,ch123data};

/* char: 0x7a 'z' */

static const unsigned char  ch122data[] = {
0xfc,0x40,0x20,0x10,0x8,0xfc,
};

static const FontGLUTBitmapCharRec ch122 = {6,6,-1,0,8,ch122data};

/* char: 0x79 'y' */

static const unsigned char  ch121data[] = {
0x78,0x84,0x4,0x74,0x8c,0x84,0x84,0x84,
};

static const FontGLUTBitmapCharRec ch121 = {6,8,-1,2,8,ch121data};

/* char: 0x78 'x' */

static const unsigned char  ch120data[] = {
0x84,0x48,0x30,0x30,0x48,0x84,
};

static const FontGLUTBitmapCharRec ch120 = {6,6,-1,0,8,ch120data};

/* char: 0x77 'w' */

static const unsigned char  ch119data[] = {
0x44,0xaa,0x92,0x92,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch119 = {7,6,0,0,8,ch119data};

/* char: 0x76 'v' */

static const unsigned char  ch118data[] = {
0x20,0x50,0x50,0x88,0x88,0x88,
};

static const FontGLUTBitmapCharRec ch118 = {5,6,-1,0,8,ch118data};

/* char: 0x75 'u' */

static const unsigned char  ch117data[] = {
0x74,0x88,0x88,0x88,0x88,0x88,
};

static const FontGLUTBitmapCharRec ch117 = {6,6,-1,0,8,ch117data};

/* char: 0x74 't' */

static const unsigned char  ch116data[] = {
0x38,0x44,0x40,0x40,0x40,0xf8,0x40,0x40,
};

static const FontGLUTBitmapCharRec ch116 = {6,8,-1,0,8,ch116data};

/* char: 0x73 's' */

static const unsigned char  ch115data[] = {
0x78,0x84,0x18,0x60,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch115 = {6,6,-1,0,8,ch115data};

/* char: 0x72 'r' */

static const unsigned char  ch114data[] = {
0x40,0x40,0x40,0x40,0x44,0xb8,
};

static const FontGLUTBitmapCharRec ch114 = {6,6,-1,0,8,ch114data};

/* char: 0x71 'q' */

static const unsigned char  ch113data[] = {
0x4,0x4,0x4,0x74,0x8c,0x84,0x8c,0x74,
};

static const FontGLUTBitmapCharRec ch113 = {6,8,-1,2,8,ch113data};

/* char: 0x70 'p' */

static const unsigned char  ch112data[] = {
0x80,0x80,0x80,0xb8,0xc4,0x84,0xc4,0xb8,
};

static const FontGLUTBitmapCharRec ch112 = {6,8,-1,2,8,ch112data};

/* char: 0x6f 'o' */

static const unsigned char  ch111data[] = {
0x78,0x84,0x84,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch111 = {6,6,-1,0,8,ch111data};

/* char: 0x6e 'n' */

static const unsigned char  ch110data[] = {
0x84,0x84,0x84,0x84,0xc4,0xb8,
};

static const FontGLUTBitmapCharRec ch110 = {6,6,-1,0,8,ch110data};

/* char: 0x6d 'm' */

static const unsigned char  ch109data[] = {
0x82,0x92,0x92,0x92,0x92,0xec,
};

static const FontGLUTBitmapCharRec ch109 = {7,6,0,0,8,ch109data};

/* char: 0x6c 'l' */

static const unsigned char  ch108data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0x60,
};

static const FontGLUTBitmapCharRec ch108 = {5,9,-1,0,8,ch108data};

/* char: 0x6b 'k' */

static const unsigned char  ch107data[] = {
0x84,0x88,0x90,0xe0,0x90,0x88,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch107 = {6,9,-1,0,8,ch107data};

/* char: 0x6a 'j' */

static const unsigned char  ch106data[] = {
0x70,0x88,0x88,0x8,0x8,0x8,0x8,0x18,0x0,0x8,
};

static const FontGLUTBitmapCharRec ch106 = {5,10,-1,2,8,ch106data};

/* char: 0x69 'i' */

static const unsigned char  ch105data[] = {
0xf8,0x20,0x20,0x20,0x20,0x60,0x0,0x20,
};

static const FontGLUTBitmapCharRec ch105 = {5,8,-1,0,8,ch105data};

/* char: 0x68 'h' */

static const unsigned char  ch104data[] = {
0x84,0x84,0x84,0x84,0xc4,0xb8,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch104 = {6,9,-1,0,8,ch104data};

/* char: 0x67 'g' */

static const unsigned char  ch103data[] = {
0x78,0x84,0x78,0x80,0x70,0x88,0x88,0x74,
};

static const FontGLUTBitmapCharRec ch103 = {6,8,-1,2,8,ch103data};

/* char: 0x66 'f' */

static const unsigned char  ch102data[] = {
0x40,0x40,0x40,0x40,0xf8,0x40,0x40,0x44,0x38,
};

static const FontGLUTBitmapCharRec ch102 = {6,9,-1,0,8,ch102data};

/* char: 0x65 'e' */

static const unsigned char  ch101data[] = {
0x78,0x84,0x80,0xfc,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch101 = {6,6,-1,0,8,ch101data};

/* char: 0x64 'd' */

static const unsigned char  ch100data[] = {
0x74,0x8c,0x84,0x84,0x8c,0x74,0x4,0x4,0x4,
};

static const FontGLUTBitmapCharRec ch100 = {6,9,-1,0,8,ch100data};

/* char: 0x63 'c' */

static const unsigned char  ch99data[] = {
0x78,0x84,0x80,0x80,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch99 = {6,6,-1,0,8,ch99data};

/* char: 0x62 'b' */

static const unsigned char  ch98data[] = {
0xb8,0xc4,0x84,0x84,0xc4,0xb8,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch98 = {6,9,-1,0,8,ch98data};

/* char: 0x61 'a' */

static const unsigned char  ch97data[] = {
0x74,0x8c,0x84,0x7c,0x4,0x78,
};

static const FontGLUTBitmapCharRec ch97 = {6,6,-1,0,8,ch97data};

/* char: 0x60 '`' */

static const unsigned char  ch96data[] = {
0x10,0x60,0xe0,
};

static const FontGLUTBitmapCharRec ch96 = {4,3,-2,-6,8,ch96data};

/* char: 0x5f '_' */

static const unsigned char  ch95data[] = {
0xfe,
};

static const FontGLUTBitmapCharRec ch95 = {7,1,0,1,8,ch95data};

/* char: 0x5e '^' */

static const unsigned char  ch94data[] = {
0x88,0x50,0x20,
};

static const FontGLUTBitmapCharRec ch94 = {5,3,-1,-6,8,ch94data};

/* char: 0x5d ']' */

static const unsigned char  ch93data[] = {
0xf0,0x10,0x10,0x10,0x10,0x10,0x10,0x10,0xf0,
};

static const FontGLUTBitmapCharRec ch93 = {4,9,-1,0,8,ch93data};

/* char: 0x5c '\' */

static const unsigned char  ch92data[] = {
0x2,0x2,0x4,0x8,0x10,0x20,0x40,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch92 = {7,9,0,0,8,ch92data};

/* char: 0x5b '[' */

static const unsigned char  ch91data[] = {
0xf0,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0xf0,
};

static const FontGLUTBitmapCharRec ch91 = {4,9,-2,0,8,ch91data};

/* char: 0x5a 'Z' */

static const unsigned char  ch90data[] = {
0xfc,0x80,0x80,0x40,0x20,0x10,0x8,0x4,0xfc,
};

static const FontGLUTBitmapCharRec ch90 = {6,9,-1,0,8,ch90data};

/* char: 0x59 'Y' */

static const unsigned char  ch89data[] = {
0x10,0x10,0x10,0x10,0x10,0x28,0x44,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch89 = {7,9,0,0,8,ch89data};

/* char: 0x58 'X' */

static const unsigned char  ch88data[] = {
0x82,0x82,0x44,0x28,0x10,0x28,0x44,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch88 = {7,9,0,0,8,ch88data};

/* char: 0x57 'W' */

static const unsigned char  ch87data[] = {
0x44,0xaa,0x92,0x92,0x92,0x82,0x82,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch87 = {7,9,0,0,8,ch87data};

/* char: 0x56 'V' */

static const unsigned char  ch86data[] = {
0x10,0x28,0x28,0x28,0x44,0x44,0x44,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch86 = {7,9,0,0,8,ch86data};

/* char: 0x55 'U' */

static const unsigned char  ch85data[] = {
0x78,0x84,0x84,0x84,0x84,0x84,0x84,0x84,0x84,
};

static const FontGLUTBitmapCharRec ch85 = {6,9,-1,0,8,ch85data};

/* char: 0x54 'T' */

static const unsigned char  ch84data[] = {
0x10,0x10,0x10,0x10,0x10,0x10,0x10,0x10,0xfe,
};

static const FontGLUTBitmapCharRec ch84 = {7,9,0,0,8,ch84data};

/* char: 0x53 'S' */

static const unsigned char  ch83data[] = {
0x78,0x84,0x4,0x4,0x78,0x80,0x80,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch83 = {6,9,-1,0,8,ch83data};

/* char: 0x52 'R' */

static const unsigned char  ch82data[] = {
0x84,0x88,0x90,0xa0,0xf8,0x84,0x84,0x84,0xf8,
};

static const FontGLUTBitmapCharRec ch82 = {6,9,-1,0,8,ch82data};

/* char: 0x51 'Q' */

static const unsigned char  ch81data[] = {
0x4,0x78,0x94,0xa4,0x84,0x84,0x84,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch81 = {6,10,-1,1,8,ch81data};

/* char: 0x50 'P' */

static const unsigned char  ch80data[] = {
0x80,0x80,0x80,0x80,0xf8,0x84,0x84,0x84,0xf8,
};

static const FontGLUTBitmapCharRec ch80 = {6,9,-1,0,8,ch80data};

/* char: 0x4f 'O' */

static const unsigned char  ch79data[] = {
0x78,0x84,0x84,0x84,0x84,0x84,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch79 = {6,9,-1,0,8,ch79data};

/* char: 0x4e 'N' */

static const unsigned char  ch78data[] = {
0x84,0x84,0x84,0x8c,0x94,0xa4,0xc4,0x84,0x84,
};

static const FontGLUTBitmapCharRec ch78 = {6,9,-1,0,8,ch78data};

/* char: 0x4d 'M' */

static const unsigned char  ch77data[] = {
0x82,0x82,0x82,0x92,0x92,0xaa,0xc6,0x82,0x82,
};

static const FontGLUTBitmapCharRec ch77 = {7,9,0,0,8,ch77data};

/* char: 0x4c 'L' */

static const unsigned char  ch76data[] = {
0xfc,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch76 = {6,9,-1,0,8,ch76data};

/* char: 0x4b 'K' */

static const unsigned char  ch75data[] = {
0x84,0x88,0x90,0xa0,0xc0,0xa0,0x90,0x88,0x84,
};

static const FontGLUTBitmapCharRec ch75 = {6,9,-1,0,8,ch75data};

/* char: 0x4a 'J' */

static const unsigned char  ch74data[] = {
0x70,0x88,0x8,0x8,0x8,0x8,0x8,0x8,0x3c,
};

static const FontGLUTBitmapCharRec ch74 = {6,9,-1,0,8,ch74data};

/* char: 0x49 'I' */

static const unsigned char  ch73data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0x20,0x20,0xf8,
};

static const FontGLUTBitmapCharRec ch73 = {5,9,-1,0,8,ch73data};

/* char: 0x48 'H' */

static const unsigned char  ch72data[] = {
0x84,0x84,0x84,0x84,0xfc,0x84,0x84,0x84,0x84,
};

static const FontGLUTBitmapCharRec ch72 = {6,9,-1,0,8,ch72data};

/* char: 0x47 'G' */

static const unsigned char  ch71data[] = {
0x74,0x8c,0x84,0x9c,0x80,0x80,0x80,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch71 = {6,9,-1,0,8,ch71data};

/* char: 0x46 'F' */

static const unsigned char  ch70data[] = {
0x80,0x80,0x80,0x80,0xf0,0x80,0x80,0x80,0xfc,
};

static const FontGLUTBitmapCharRec ch70 = {6,9,-1,0,8,ch70data};

/* char: 0x45 'E' */

static const unsigned char  ch69data[] = {
0xfc,0x80,0x80,0x80,0xf0,0x80,0x80,0x80,0xfc,
};

static const FontGLUTBitmapCharRec ch69 = {6,9,-1,0,8,ch69data};

/* char: 0x44 'D' */

static const unsigned char  ch68data[] = {
0xfc,0x42,0x42,0x42,0x42,0x42,0x42,0x42,0xfc,
};

static const FontGLUTBitmapCharRec ch68 = {7,9,0,0,8,ch68data};

/* char: 0x43 'C' */

static const unsigned char  ch67data[] = {
0x78,0x84,0x80,0x80,0x80,0x80,0x80,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch67 = {6,9,-1,0,8,ch67data};

/* char: 0x42 'B' */

static const unsigned char  ch66data[] = {
0xfc,0x42,0x42,0x42,0x7c,0x42,0x42,0x42,0xfc,
};

static const FontGLUTBitmapCharRec ch66 = {7,9,0,0,8,ch66data};

/* char: 0x41 'A' */

static const unsigned char  ch65data[] = {
0x84,0x84,0x84,0xfc,0x84,0x84,0x84,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch65 = {6,9,-1,0,8,ch65data};

/* char: 0x40 '@' */

static const unsigned char  ch64data[] = {
0x78,0x80,0x94,0xac,0xa4,0x9c,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch64 = {6,9,-1,0,8,ch64data};

/* char: 0x3f '?' */

static const unsigned char  ch63data[] = {
0x10,0x0,0x10,0x10,0x8,0x4,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch63 = {6,9,-1,0,8,ch63data};

/* char: 0x3e '>' */

static const unsigned char  ch62data[] = {
0x80,0x40,0x20,0x10,0x8,0x10,0x20,0x40,0x80,
};

static const FontGLUTBitmapCharRec ch62 = {5,9,-1,0,8,ch62data};

/* char: 0x3d '=' */

static const unsigned char  ch61data[] = {
0xfc,0x0,0x0,0xfc,
};

static const FontGLUTBitmapCharRec ch61 = {6,4,-1,-2,8,ch61data};

/* char: 0x3c '<' */

static const unsigned char  ch60data[] = {
0x8,0x10,0x20,0x40,0x80,0x40,0x20,0x10,0x8,
};

static const FontGLUTBitmapCharRec ch60 = {5,9,-2,0,8,ch60data};

/* char: 0x3b ';' */

static const unsigned char  ch59data[] = {
0x80,0x60,0x70,0x0,0x0,0x20,0x70,0x20,
};

static const FontGLUTBitmapCharRec ch59 = {4,8,-1,1,8,ch59data};

/* char: 0x3a ':' */

static const unsigned char  ch58data[] = {
0x40,0xe0,0x40,0x0,0x0,0x40,0xe0,0x40,
};

static const FontGLUTBitmapCharRec ch58 = {3,8,-2,1,8,ch58data};

/* char: 0x39 '9' */

static const unsigned char  ch57data[] = {
0x70,0x8,0x4,0x4,0x74,0x8c,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch57 = {6,9,-1,0,8,ch57data};

/* char: 0x38 '8' */

static const unsigned char  ch56data[] = {
0x78,0x84,0x84,0x84,0x78,0x84,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch56 = {6,9,-1,0,8,ch56data};

/* char: 0x37 '7' */

static const unsigned char  ch55data[] = {
0x40,0x40,0x20,0x20,0x10,0x10,0x8,0x4,0xfc,
};

static const FontGLUTBitmapCharRec ch55 = {6,9,-1,0,8,ch55data};

/* char: 0x36 '6' */

static const unsigned char  ch54data[] = {
0x78,0x84,0x84,0xc4,0xb8,0x80,0x80,0x40,0x38,
};

static const FontGLUTBitmapCharRec ch54 = {6,9,-1,0,8,ch54data};

/* char: 0x35 '5' */

static const unsigned char  ch53data[] = {
0x78,0x84,0x4,0x4,0xc4,0xb8,0x80,0x80,0xfc,
};

static const FontGLUTBitmapCharRec ch53 = {6,9,-1,0,8,ch53data};

/* char: 0x34 '4' */

static const unsigned char  ch52data[] = {
0x8,0x8,0xfc,0x88,0x88,0x48,0x28,0x18,0x8,
};

static const FontGLUTBitmapCharRec ch52 = {6,9,-1,0,8,ch52data};

/* char: 0x33 '3' */

static const unsigned char  ch51data[] = {
0x78,0x84,0x4,0x4,0x38,0x10,0x8,0x4,0xfc,
};

static const FontGLUTBitmapCharRec ch51 = {6,9,-1,0,8,ch51data};

/* char: 0x32 '2' */

static const unsigned char  ch50data[] = {
0xfc,0x80,0x40,0x30,0x8,0x4,0x84,0x84,0x78,
};

static const FontGLUTBitmapCharRec ch50 = {6,9,-1,0,8,ch50data};

/* char: 0x31 '1' */

static const unsigned char  ch49data[] = {
0xf8,0x20,0x20,0x20,0x20,0x20,0xa0,0x60,0x20,
};

static const FontGLUTBitmapCharRec ch49 = {5,9,-1,0,8,ch49data};

/* char: 0x30 '0' */

static const unsigned char  ch48data[] = {
0x30,0x48,0x84,0x84,0x84,0x84,0x84,0x48,0x30,
};

static const FontGLUTBitmapCharRec ch48 = {6,9,-1,0,8,ch48data};

/* char: 0x2f '/' */

static const unsigned char  ch47data[] = {
0x80,0x80,0x40,0x20,0x10,0x8,0x4,0x2,0x2,
};

static const FontGLUTBitmapCharRec ch47 = {7,9,0,0,8,ch47data};

/* char: 0x2e '.' */

static const unsigned char  ch46data[] = {
0x40,0xe0,0x40,
};

static const FontGLUTBitmapCharRec ch46 = {3,3,-2,1,8,ch46data};

/* char: 0x2d '-' */

static const unsigned char  ch45data[] = {
0xfc,
};

static const FontGLUTBitmapCharRec ch45 = {6,1,-1,-4,8,ch45data};

/* char: 0x2c ',' */

static const unsigned char  ch44data[] = {
0x80,0x60,0x70,
};

static const FontGLUTBitmapCharRec ch44 = {4,3,-1,1,8,ch44data};

/* char: 0x2b '+' */

static const unsigned char  ch43data[] = {
0x20,0x20,0xf8,0x20,0x20,
};

static const FontGLUTBitmapCharRec ch43 = {5,5,-1,-2,8,ch43data};

/* char: 0x2a '*' */

static const unsigned char  ch42data[] = {
0x48,0x30,0xfc,0x30,0x48,
};

static const FontGLUTBitmapCharRec ch42 = {6,5,-1,-2,8,ch42data};

/* char: 0x29 ')' */

static const unsigned char  ch41data[] = {
0x80,0x40,0x40,0x20,0x20,0x20,0x40,0x40,0x80,
};

static const FontGLUTBitmapCharRec ch41 = {3,9,-2,0,8,ch41data};

/* char: 0x28 '(' */

static const unsigned char  ch40data[] = {
0x20,0x40,0x40,0x80,0x80,0x80,0x40,0x40,0x20,
};

static const FontGLUTBitmapCharRec ch40 = {3,9,-3,0,8,ch40data};

/* char: 0x27 ''' */

static const unsigned char  ch39data[] = {
0x80,0x60,0x70,
};

static const FontGLUTBitmapCharRec ch39 = {4,3,-1,-6,8,ch39data};

/* char: 0x26 '&' */

static const unsigned char  ch38data[] = {
0x74,0x88,0x94,0x60,0x90,0x90,0x60,
};

static const FontGLUTBitmapCharRec ch38 = {6,7,-1,0,8,ch38data};

/* char: 0x25 '%' */

static const unsigned char  ch37data[] = {
0x88,0x54,0x48,0x20,0x10,0x10,0x48,0xa4,0x44,
};

static const FontGLUTBitmapCharRec ch37 = {6,9,-1,0,8,ch37data};

/* char: 0x24 '$' */

static const unsigned char  ch36data[] = {
0x20,0xf0,0x28,0x70,0xa0,0x78,0x20,
};

static const FontGLUTBitmapCharRec ch36 = {5,7,-1,-1,8,ch36data};

/* char: 0x23 '#' */

static const unsigned char  ch35data[] = {
0x48,0x48,0xfc,0x48,0xfc,0x48,0x48,
};

static const FontGLUTBitmapCharRec ch35 = {6,7,-1,-1,8,ch35data};

/* char: 0x22 '"' */

static const unsigned char  ch34data[] = {
0x90,0x90,0x90,
};

static const FontGLUTBitmapCharRec ch34 = {4,3,-2,-6,8,ch34data};

/* char: 0x21 '!' */

static const unsigned char  ch33data[] = {
0x80,0x0,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch33 = {1,9,-3,0,8,ch33data};

/* char: 0x1f */

static const unsigned char  ch31data[] = {
0x80,
};

static const FontGLUTBitmapCharRec ch31 = {1,1,-3,-3,8,ch31data};

/* char: 0x1e */

static const unsigned char  ch30data[] = {
0xdc,0x62,0x20,0x20,0x20,0x70,0x20,0x22,0x1c,
};

static const FontGLUTBitmapCharRec ch30 = {7,9,0,0,8,ch30data};

/* char: 0x1d */

static const unsigned char  ch29data[] = {
0x80,0x40,0xfe,0x10,0xfe,0x4,0x2,
};

static const FontGLUTBitmapCharRec ch29 = {7,7,0,0,8,ch29data};

/* char: 0x1c */

static const unsigned char  ch28data[] = {
0x88,0x48,0x48,0x48,0x48,0xfc,
};

static const FontGLUTBitmapCharRec ch28 = {6,6,-1,0,8,ch28data};

/* char: 0x1b */

static const unsigned char  ch27data[] = {
0xfe,0x80,0x20,0x8,0x2,0x8,0x20,0x80,
};

static const FontGLUTBitmapCharRec ch27 = {7,8,0,0,8,ch27data};

/* char: 0x1a */

static const unsigned char  ch26data[] = {
0xfe,0x2,0x8,0x20,0x80,0x20,0x8,0x2,
};

static const FontGLUTBitmapCharRec ch26 = {7,8,0,0,8,ch26data};

/* char: 0x19 */

static const unsigned char  ch25data[] = {
0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch25 = {1,13,-3,2,8,ch25data};

/* char: 0x18 */

static const unsigned char  ch24data[] = {
0x10,0x10,0x10,0x10,0x10,0xff,
};

static const FontGLUTBitmapCharRec ch24 = {8,6,0,2,8,ch24data};

/* char: 0x17 */

static const unsigned char  ch23data[] = {
0xff,0x10,0x10,0x10,0x10,0x10,0x10,0x10,
};

static const FontGLUTBitmapCharRec ch23 = {8,8,0,-3,8,ch23data};

/* char: 0x16 */

static const unsigned char  ch22data[] = {
0x10,0x10,0x10,0x10,0x10,0xf0,0x10,0x10,0x10,0x10,0x10,0x10,0x10,
};

static const FontGLUTBitmapCharRec ch22 = {4,13,0,2,8,ch22data};

/* char: 0x15 */

static const unsigned char  ch21data[] = {
0x80,0x80,0x80,0x80,0x80,0xf8,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch21 = {5,13,-3,2,8,ch21data};

/* char: 0x14 */

static const unsigned char  ch20data[] = {
0xff,
};

static const FontGLUTBitmapCharRec ch20 = {8,1,0,1,8,ch20data};

/* char: 0x13 */

static const unsigned char  ch19data[] = {
0xff,
};

static const FontGLUTBitmapCharRec ch19 = {8,1,0,-1,8,ch19data};

/* char: 0x12 */

static const unsigned char  ch18data[] = {
0xff,
};

static const FontGLUTBitmapCharRec ch18 = {8,1,0,-3,8,ch18data};

/* char: 0x11 */

static const unsigned char  ch17data[] = {
0xff,
};

static const FontGLUTBitmapCharRec ch17 = {8,1,0,-5,8,ch17data};

/* char: 0x10 */

static const unsigned char  ch16data[] = {
0xff,
};

static const FontGLUTBitmapCharRec ch16 = {8,1,0,-7,8,ch16data};

/* char: 0xf */

static const unsigned char  ch15data[] = {
0x10,0x10,0x10,0x10,0x10,0xff,0x10,0x10,0x10,0x10,0x10,0x10,0x10,
};

static const FontGLUTBitmapCharRec ch15 = {8,13,0,2,8,ch15data};

/* char: 0xe */

static const unsigned char  ch14data[] = {
0xf8,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch14 = {5,8,-3,-3,8,ch14data};

/* char: 0xd */

static const unsigned char  ch13data[] = {
0x80,0x80,0x80,0x80,0x80,0xf8,
};

static const FontGLUTBitmapCharRec ch13 = {5,6,-3,2,8,ch13data};

/* char: 0xc */

static const unsigned char  ch12data[] = {
0x10,0x10,0x10,0x10,0x10,0xf0,
};

static const FontGLUTBitmapCharRec ch12 = {4,6,0,2,8,ch12data};

/* char: 0xb */

static const unsigned char  ch11data[] = {
0xf0,0x10,0x10,0x10,0x10,0x10,0x10,0x10,
};

static const FontGLUTBitmapCharRec ch11 = {4,8,0,-3,8,ch11data};

/* char: 0xa */

static const unsigned char  ch10data[] = {
0x8,0x8,0x8,0x8,0x3e,0x20,0x50,0x88,0x88,
};

static const FontGLUTBitmapCharRec ch10 = {7,9,0,2,8,ch10data};

/* char: 0x9 */

static const unsigned char  ch9data[] = {
0x3e,0x20,0x20,0x20,0x88,0x98,0xa8,0xc8,0x88,
};

static const FontGLUTBitmapCharRec ch9 = {7,9,0,2,8,ch9data};

/* char: 0x8 */

static const unsigned char  ch8data[] = {
0xfe,0x10,0x10,0xfe,0x10,0x10,
};

static const FontGLUTBitmapCharRec ch8 = {7,6,0,0,8,ch8data};

/* char: 0x7 */

static const unsigned char  ch7data[] = {
0x70,0x88,0x88,0x70,
};

static const FontGLUTBitmapCharRec ch7 = {5,4,-1,-5,8,ch7data};

/* char: 0x6 */

static const unsigned char  ch6data[] = {
0x20,0x20,0x3c,0x20,0x3e,0xf8,0x80,0x80,0x80,
};

static const FontGLUTBitmapCharRec ch6 = {7,9,0,2,8,ch6data};

/* char: 0x5 */

static const unsigned char  ch5data[] = {
0x22,0x22,0x3c,0x22,0x3c,0x78,0x80,0x80,0x78,
};

static const FontGLUTBitmapCharRec ch5 = {7,9,0,2,8,ch5data};

/* char: 0x4 */

static const unsigned char  ch4data[] = {
0x10,0x10,0x1c,0x10,0x9e,0x80,0xe0,0x80,0xf0,
};

static const FontGLUTBitmapCharRec ch4 = {7,9,0,2,8,ch4data};

/* char: 0x3 */

static const unsigned char  ch3data[] = {
0x8,0x8,0x8,0x3e,0x88,0x88,0xf8,0x88,0x88,
};

static const FontGLUTBitmapCharRec ch3 = {7,9,0,2,8,ch3data};

/* char: 0x2 */

static const unsigned char  ch2data[] = {
0x55,0xaa,0x55,0xaa,0x55,0xaa,0x55,0xaa,0x55,0xaa,0x55,0xaa,
};

static const FontGLUTBitmapCharRec ch2 = {8,12,0,2,8,ch2data};

/* char: 0x1 */

static const unsigned char  ch1data[] = {
0x10,0x38,0x7c,0xfe,0x7c,0x38,0x10,
};

static const FontGLUTBitmapCharRec ch1 = {7,7,0,-1,8,ch1data};

static const FontGLUTBitmapCharRec * const chars[] = {
&ch0,
&ch1,
&ch2,
&ch3,
&ch4,
&ch5,
&ch6,
&ch7,
&ch8,
&ch9,
&ch10,
&ch11,
&ch12,
&ch13,
&ch14,
&ch15,
&ch16,
&ch17,
&ch18,
&ch19,
&ch20,
&ch21,
&ch22,
&ch23,
&ch24,
&ch25,
&ch26,
&ch27,
&ch28,
&ch29,
&ch30,
&ch31,
&ch32,
&ch33,
&ch34,
&ch35,
&ch36,
&ch37,
&ch38,
&ch39,
&ch40,
&ch41,
&ch42,
&ch43,
&ch44,
&ch45,
&ch46,
&ch47,
&ch48,
&ch49,
&ch50,
&ch51,
&ch52,
&ch53,
&ch54,
&ch55,
&ch56,
&ch57,
&ch58,
&ch59,
&ch60,
&ch61,
&ch62,
&ch63,
&ch64,
&ch65,
&ch66,
&ch67,
&ch68,
&ch69,
&ch70,
&ch71,
&ch72,
&ch73,
&ch74,
&ch75,
&ch76,
&ch77,
&ch78,
&ch79,
&ch80,
&ch81,
&ch82,
&ch83,
&ch84,
&ch85,
&ch86,
&ch87,
&ch88,
&ch89,
&ch90,
&ch91,
&ch92,
&ch93,
&ch94,
&ch95,
&ch96,
&ch97,
&ch98,
&ch99,
&ch100,
&ch101,
&ch102,
&ch103,
&ch104,
&ch105,
&ch106,
&ch107,
&ch108,
&ch109,
&ch110,
&ch111,
&ch112,
&ch113,
&ch114,
&ch115,
&ch116,
&ch117,
&ch118,
&ch119,
&ch120,
&ch121,
&ch122,
&ch123,
&ch124,
&ch125,
&ch126,
&ch127,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
&ch160,
&ch161,
&ch162,
&ch163,
&ch164,
&ch165,
&ch166,
&ch167,
&ch168,
&ch169,
&ch170,
&ch171,
&ch172,
&ch173,
&ch174,
&ch175,
&ch176,
&ch177,
&ch178,
&ch179,
&ch180,
&ch181,
&ch182,
&ch183,
&ch184,
&ch185,
&ch186,
&ch187,
&ch188,
&ch189,
&ch190,
&ch191,
&ch192,
&ch193,
&ch194,
&ch195,
&ch196,
&ch197,
&ch198,
&ch199,
&ch200,
&ch201,
&ch202,
&ch203,
&ch204,
&ch205,
&ch206,
&ch207,
&ch208,
&ch209,
&ch210,
&ch211,
&ch212,
&ch213,
&ch214,
&ch215,
&ch216,
&ch217,
&ch218,
&ch219,
&ch220,
&ch221,
&ch222,
&ch223,
&ch224,
&ch225,
&ch226,
&ch227,
&ch228,
&ch229,
&ch230,
&ch231,
&ch232,
&ch233,
&ch234,
&ch235,
&ch236,
&ch237,
&ch238,
&ch239,
&ch240,
&ch241,
&ch242,
&ch243,
&ch244,
&ch245,
&ch246,
&ch247,
&ch248,
&ch249,
&ch250,
&ch251,
&ch252,
&ch253,
&ch254,
&ch255,
};

FontGLUTBitmapFontRec FontGLUTBitmap8By13 = {
"-misc-fixed-medium-r-normal--13-120-75-75-C-80-iso8859-1",
256,
0,
chars
};

