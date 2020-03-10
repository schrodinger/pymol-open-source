#pragma once

/**
 * Enums for CGOEnable
 */
enum {
  GL_DEFAULT_SHADER_WITH_SETTINGS       = 0xffe0,
  GL_SPHERE_SHADER                      = 0xffe1,
  GL_CYLINDER_SHADER                    = 0xffe2,
  GL_TWO_SIDED_LIGHTING                 = 0xffe3,
  GL_MESH_LIGHTING                      = 0xffe4,
  GL_DOT_LIGHTING                       = 0xffe5,
  GL_LABEL_FLOAT_TEXT                   = 0xffe6,
  GL_DASH_TRANSPARENCY_DEPTH_TEST       = 0xffe7,
  GL_BACK_FACE_CULLING                  = 0xffe8,
  GL_DEPTH_TEST_IF_FLOATING             = 0xffe9,
  GL_OIT_COPY_SHADER                    = 0xffea,
  GL_SURFACE_SHADER                     = 0xffeb,
  GL_LINE_SHADER                        = 0xffec,

  CGO_GL_LIGHTING                       = 0xffef,

  GL_SCREEN_SHADER                      = 0xfff1,
  GL_RAMP_SHADER                        = 0xfff2,
  GL_CONNECTOR_SHADER                   = 0xfff3,
  GL_FXAA_SHADER                        = 0xfff4,
  GL_SMAA1_SHADER                       = 0xfff5,
  GL_SMAA2_SHADER                       = 0xfff6,
  GL_SMAA3_SHADER                       = 0xfff7,
  GL_TRILINES_SHADER                    = 0xfff8,
  GL_OIT_SHADER                         = 0xfff9,
  GL_LABEL_SHADER                       = 0xfffa,
  GL_BACKGROUND_SHADER                  = 0xfffb,

  GL_DEFAULT_SHADER                     = 0xfffd,
  GL_SHADER_LIGHTING                    = 0xfffe,
};
