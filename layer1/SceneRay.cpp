
#include"Scene.h"
#include"SceneRay.h"
#include"SceneDef.h"
#include"Util.h"
#include"ShaderMgr.h"
#include"Matrix.h"
#include"PyMOL.h"
#include"ListMacros.h"
#include"Color.h"
#include"P.h"

static double accumTiming = 0.0;

/* EXPERIMENTAL VOLUME RAYTRACING DATA */
static std::shared_ptr<pymol::Image> rayVolumeImage;
extern float *rayDepthPixels;
extern int rayVolume, rayWidth, rayHeight;


static void SceneRaySetRayView(PyMOLGlobals * G, CScene *I, int stereo_hand,
    float *rayView, float *angle, float shift)
{
  /* start afresh, looking in the negative Z direction (0,0,-1) from (0,0,0) */
  identity44f(rayView);
  
  if(stereo_hand) {
    /* stereo */
    
    float stAng, stShift;
    stAng = SettingGetGlobal_f(G, cSetting_stereo_angle);
    stShift = SettingGetGlobal_f(G, cSetting_stereo_shift);
    /* right hand */
    stShift = (float) (stShift * fabs(I->m_view.m_pos[2]) / 100.0);
    stAng = (float) (stAng * atan(stShift / fabs(I->m_view.m_pos[2])) * 90.0 / cPI);
    if(stereo_hand == 2) {  /* left hand */
      stAng = -stAng;
      stShift = -stShift;
    }
    *angle = stAng;
    {
      float temp[16];
      identity44f(temp);
      MatrixRotateC44f(temp, (float) (-PI * stAng / 180), 0.0F, 1.0F, 0.0F);        /* y-axis rotation */
      MatrixMultiplyC44f(temp, rayView);
    }
    /* move the camera to the location we are looking at */
    MatrixTranslateC44f(rayView, I->m_view.m_pos[0], I->m_view.m_pos[1], I->m_view.m_pos[2]);
    MatrixTranslateC44f(rayView, stShift, 0.0, 0.0);
    MatrixMultiplyC44f(I->m_view.m_rotMatrix, rayView);
  } else {                  /* not stereo mode */
    /* move the camera to the location we are looking at */
    MatrixTranslateC44f(rayView, I->m_view.m_pos[0], I->m_view.m_pos[1], I->m_view.m_pos[2]);
    if(shift) {
      MatrixTranslateC44f(rayView, shift, 0.0F, 0.0F);
    }
    /* move the camera so that we can see the origin 
     * NOTE, vector is given in the coordinates of the world's motion
     * relative to the camera */
    /* 4. rotate about the origin (the the center of rotation) */
    if(*angle) {
      float temp[16];
      identity44f(temp);
      MatrixRotateC44f(temp, (float) (-PI * *angle / 180), 0.0F, 1.0F, 0.0F);
      MatrixMultiplyC44f(I->m_view.m_rotMatrix, temp);
      MatrixMultiplyC44f(temp, rayView);
    } else {
      MatrixMultiplyC44f(I->m_view.m_rotMatrix, rayView);
    }
  }
  /* 5. move the origin to the center of rotation */
  MatrixTranslateC44f(rayView, -I->m_view.m_origin[0], -I->m_view.m_origin[1], -I->m_view.m_origin[2]);
  
  if(Feedback(G, FB_Scene, FB_Debugging)) {
    fprintf(stderr, "SceneRay: %8.3f %8.3f %8.3f\n", I->m_view.m_pos[0], I->m_view.m_pos[1], I->m_view.m_pos[2]);
    fprintf(stderr, "SceneRay: %8.3f %8.3f %8.3f\n",
	    I->m_view.m_origin[0], I->m_view.m_origin[1], I->m_view.m_origin[2]);
    fprintf(stderr, "SceneRay: %8.3f %8.3f %8.3f\n",
	    I->m_view.m_rotMatrix[0], I->m_view.m_rotMatrix[1], I->m_view.m_rotMatrix[2]);
  }
}

bool SceneRay(PyMOLGlobals * G,
              int ray_width, int ray_height, int mode,
              char **headerVLA_ptr,
              char **charVLA_ptr, float angle,
              float shift, int quiet, G3dPrimitive ** g3d, int show_timing, int antialias)
{
#ifdef _PYMOL_NO_RAY
  FeedbackAdd(G, "" _PYMOL_NO_RAY);
  return false;
#else

  CScene *I = G->Scene;
  CRay *ray = NULL;
  float height, width;
  float aspRat;
  float rayView[16];
  double timing;
  char *charVLA = NULL;
  char *headerVLA = NULL;
  float fov;
  int stereo_hand = 0;
  int grid_mode = SettingGetGlobal_i(G, cSetting_grid_mode);
  std::shared_ptr<pymol::Image> stereo_image;
  OrthoLineType prefix = "";
  int ortho = SettingGetGlobal_i(G, cSetting_ray_orthoscopic);
  int last_grid_active = I->grid.active;
  int grid_size = 0;

  if(SettingGetGlobal_i(G, cSetting_defer_builds_mode) == 5)
    SceneUpdate(G, true);

  if(ortho < 0)
    ortho = SettingGetGlobal_b(G, cSetting_ortho);

  if(mode != 0)
    grid_mode = 0;              /* only allow grid mode with PyMOL renderer */

  SceneUpdateAnimation(G);

  if(mode == 0)
    SceneInvalidateCopy(G, true);

  if(antialias < 0) {
    antialias = SettingGetGlobal_i(G, cSetting_antialias);
  }
  if(ray_width < 0)
    ray_width = 0;
  if(ray_height < 0)
    ray_height = 0;
  if((!ray_width) || (!ray_height)) {
    if(ray_width && (!ray_height)) {
      ray_height = (ray_width * I->Height) / I->Width;
    } else if(ray_height && (!ray_width)) {
      ray_width = (ray_height * I->Width) / I->Height;
    } else {
      ray_width = I->Width;
      ray_height = I->Height;
    }
  }

  fov = SettingGetGlobal_f(G, cSetting_field_of_view);

  timing = UtilGetSeconds(G);   /* start timing the process */

  SceneUpdate(G, false);

  switch (I->StereoMode) {
  case cStereo_quadbuffer:
  case cStereo_openvr:
    stereo_hand = 2;
    break;
  case cStereo_crosseye:
  case cStereo_walleye:
    ray_width = ray_width / 2;
    stereo_hand = 2;
    break;
  case cStereo_geowall:
  case cStereo_sidebyside:
    stereo_hand = 2;
    break;
  case cStereo_stencil_by_row:
  case cStereo_stencil_by_column:
  case cStereo_stencil_checkerboard:
  case cStereo_stencil_custom:
  case cStereo_anaglyph:
    stereo_hand = 2;
    break;
  }

  aspRat = ((float) ray_width) / ((float) ray_height);

  if(grid_mode) {
    grid_size = SceneGetGridSize(G, grid_mode);
    GridUpdate(&I->grid, aspRat, grid_mode, grid_size);
    if(I->grid.active)
      aspRat *= I->grid.asp_adjust;
  }
  if (last_grid_active != I->grid.active || grid_size != I->last_grid_size){
    //    ExecutiveInvalidateRep(G, cKeywordAll, cRepLabel, cRepInvAll);
    G->ShaderMgr->ResetUniformSet();    
  }
  I->last_grid_size = grid_size;
  while(1) {
    int slot;
    int tot_width = ray_width;
    int tot_height = ray_height;
    int ray_x = 0, ray_y = 0;

    if(I->grid.active)
      GridGetRayViewport(&I->grid, ray_width, ray_height);

    for(slot = 0; slot <= I->grid.last_slot; slot++) {

      if(I->grid.active) {
        GridSetRayViewport(&I->grid, slot, &ray_x, &ray_y, &ray_width, &ray_height);
        OrthoBusySlow(G, slot, I->grid.last_slot);
      }

      ray = RayNew(G, antialias);
      if(!ray)
        break;

      SceneRaySetRayView(G, I, stereo_hand, rayView, &angle, shift);

      /* define the viewing volume */

      height = (float) (fabs(I->m_view.m_pos[2]) * tan((fov / 2.0) * cPI / 180.0));
      width = height * aspRat;
      PyMOL_SetBusy(G->PyMOL, true);
      OrthoBusyFast(G, 0, 20);

      {
        float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);

        if(pixel_scale_value < 0)
          pixel_scale_value = 1.0F;

        pixel_scale_value *= ((float) tot_height) / I->Height;

        if(ortho) {
          const float _1 = 1.0F;
          RayPrepare(ray, -width, width, -height, height, I->m_view.m_clipSafe.m_front,
                     I->m_view.m_clipSafe.m_back, fov,  I->m_view.m_pos, rayView, I->m_view.m_rotMatrix,
                     aspRat, ray_width, ray_height, 
                     pixel_scale_value, ortho, _1, _1,      
                     ((float) ray_height) / I->Height);
        } else {
          float back_ratio;
          float back_height;
          float back_width;
          float pos;
          pos = I->m_view.m_pos[2];

          if((-pos) < I->m_view.m_clipSafe.m_front) {
            pos = -I->m_view.m_clipSafe.m_front;
          }

          back_ratio = -I->m_view.m_clipSafe.m_back / pos;
          back_height = back_ratio * height;
          back_width = aspRat * back_height;
          RayPrepare(ray,
                     -back_width, back_width,
                     -back_height, back_height,
                     I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back,
                     fov, I->m_view.m_pos,
                     rayView, I->m_view.m_rotMatrix, aspRat,
                     ray_width, ray_height,
                     pixel_scale_value, ortho,
                     height / back_height,
                     I->m_view.m_clipSafe.m_front / I->m_view.m_clipSafe.m_back, ((float) ray_height) / I->Height);
        }
      }
      {
        int *slot_vla = I->SlotVLA;
        int state = SceneGetState(G);
        RenderInfo info;
        info.ray = ray;
        info.ortho = ortho;
        info.vertex_scale = SceneGetScreenVertexScale(G, NULL);
	info.use_shaders = SettingGetGlobal_b(G, cSetting_use_shaders);

        if(SettingGetGlobal_b(G, cSetting_dynamic_width)) {
          info.dynamic_width = true;
          info.dynamic_width_factor =
            SettingGetGlobal_f(G, cSetting_dynamic_width_factor);
          info.dynamic_width_min = SettingGetGlobal_f(G, cSetting_dynamic_width_min);
          info.dynamic_width_max = SettingGetGlobal_f(G, cSetting_dynamic_width_max);
        }

        for (auto* obj : I->Obj) {
          // ObjectGroup used to have fRender = NULL
          if (obj->type != cObjectGroup) {
            if(SceneGetDrawFlag(&I->grid, slot_vla, obj->grid_slot)) {
              float color[3];
              ColorGetEncoded(G, obj->Color, color);
              RaySetContext(ray, obj->getRenderContext());
              ray->color3fv(color);

              auto icx = SettingGetWD<int>(
                  obj->Setting.get(), cSetting_ray_interior_color, cColorDefault);

              if (icx == cColorDefault) {
                ray->interiorColor3fv(color, true);
              } else if (icx == cColorObject) {
                ray->interiorColor3fv(color, false);
              } else {
                float icolor[3];
                ColorGetEncoded(G, icx, icolor);
                ray->interiorColor3fv(icolor, false);
              }

              if((!I->grid.active) || (I->grid.mode < 2)) {
                info.state = ObjectGetCurrentState(obj, false);
                obj->render(&info);
              } else if(I->grid.slot) {
                if (I->grid.mode == 2) {
                  if((info.state = state + I->grid.slot - 1) >= 0)
                    obj->render(&info);
                } else if (I->grid.mode == 3) {
                  info.state = I->grid.slot - obj->grid_slot - 1;
                  if (info.state >= 0 && info.state < obj->getNFrame())
                    obj->render(&info);
                }
              }
            }
          }
        }
      }

      OrthoBusyFast(G, 1, 20);

      if(mode != 2) {           /* don't show pixel count for tests */
        if(!quiet) {
          PRINTFB(G, FB_Ray, FB_Blather)
            " Ray: tracing %dx%d = %d rays against %d primitives.\n", ray_width,
            ray_height, ray_width * ray_height, RayGetNPrimitives(ray)
            ENDFB(G);
        }
      }
      switch (mode) {
      case 0:                  /* mode 0 is built-in */
        {
          auto image = pymol::make_unique<pymol::Image>(ray_width, ray_height);
          std::uint32_t background;

          RayRender(ray, image->pixels(), timing, angle, antialias, &background);

          /*    RayRenderColorTable(ray,ray_width,ray_height,buffer); */
          if(!I->grid.active) {
            I->Image = std::move(image);
          } else {
            if(!I->Image) {     /* alloc on first pass */
              I->Image = pymol::make_unique<pymol::Image>(tot_width, tot_height);
              if(I->Image) {
                unsigned int tot_size = tot_width * tot_height;
                {               /* fill with background color */
                  unsigned int *ptr = I->Image->pixels();
                  for(size_t i = 0; i < tot_size; ++i) {
                    *(ptr++) = background;
                  }
                }
              }
            }
            /* merge in the latest rendering */
            if(I->Image && I->Image->bits()) {
              int i, j;
              unsigned int *src = image->pixels();
              unsigned int *dst = I->Image->pixels();

              dst += (ray_x + ray_y * tot_width);

              for(i = 0; i < ray_height; i++) {
                for(j = 0; j < ray_width; j++) {
                  if(*src != background)
                    *(dst) = *(src);
                  dst++;
                  src++;
                }
                dst += (tot_width - ray_width);
              }
            }
          }
          I->DirtyFlag = false;
          I->CopyType = true;
          I->CopyForced = true;

          if (SettingGet<bool>(G, cSetting_ray_volume) && !I->Image->empty()) {
            rayVolumeImage = I->Image;
          } else {
            rayVolumeImage = nullptr;
          }
        }
        break;

      case 1:                  /* mode 1 is povray */
        charVLA = VLACalloc(char, 100000);
        headerVLA = VLACalloc(char, 2000);
        RayRenderPOV(ray, ray_width, ray_height, &headerVLA, &charVLA,
                     I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov, angle, antialias);
        if(!(charVLA_ptr && headerVLA_ptr)) {   /* immediate mode */
          strcpy(prefix, SettingGet_s(G, NULL, NULL, cSetting_batch_prefix));
#ifndef _PYMOL_NOPY
          if(PPovrayRender(G, headerVLA, charVLA, prefix, ray_width,
                           ray_height, antialias)) {
            strcat(prefix, ".png");
            SceneLoadPNG(G, prefix, false, 0, false);
            I->DirtyFlag = false;
          }
#endif
          VLAFreeP(charVLA);
          VLAFreeP(headerVLA);
        } else {                /* get_povray mode */
          *charVLA_ptr = charVLA;
          *headerVLA_ptr = headerVLA;
        }
        break;
      case 2:                  /* mode 2 is for testing of geometries */
        RayRenderTest(ray, ray_width, ray_height, I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov);
        break;
      case 3:                  /* mode 3 is for Jmol */
        {
          G3dPrimitive *jp =
            RayRenderG3d(ray, ray_width, ray_height, I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov,
                         quiet);
          if(0) {
            int cnt = VLAGetSize(jp);
            int a;
            for(a = 0; a < cnt; a++) {
              switch (jp[a].op) {
              case 1:
                printf("g3d.fillSphereCentered(gray,%d,%d,%d,%d);\n", jp[a].r, jp[a].x1,
                       jp[a].y1, jp[a].z1);
                break;
              case 2:
                printf("triangle(%d,%d,%d,%d,%d,%d,%d,%d,%d);\n",
                       jp[a].x1, jp[a].y1, jp[a].z1,
                       jp[a].x2, jp[a].y2, jp[a].z2, jp[a].x3, jp[a].y3, jp[a].z3);
                break;
              case 3:
                printf("g3d.fillCylinder(gray,gray,(byte)3,%d,%d,%d,%d,%d,%d,%d);\n",
                       jp[a].r,
                       jp[a].x1, jp[a].y1, jp[a].z1, jp[a].x2, jp[a].y2, jp[a].z2);
                break;
              }
            }
          }
          if(g3d) {
            *g3d = jp;
          } else {
            VLAFreeP(jp);
          }
        }
        break;
      case 4:                  /* VRML2 */
        {
          char *vla = VLACalloc(char, 100000);
          RayRenderVRML2(ray, ray_width, ray_height, &vla,
                         I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov, angle, I->m_view.m_pos[2]);
          *charVLA_ptr = vla;
        }
        break;
      case 5:                  /* mode 5 is OBJ MTL */
        {
          char *objVLA = VLACalloc(char, 100000);
          char *mtlVLA = VLACalloc(char, 1000);
          RayRenderObjMtl(ray, ray_width, ray_height, &objVLA, &mtlVLA,
                          I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov, angle, I->m_view.m_pos[2]);
          *headerVLA_ptr = objVLA;
          *charVLA_ptr = mtlVLA;
        }
        break;
      case 6:                  /* VRML1 -- more compatible with tools like blender */
        {
          char *vla = VLACalloc(char, 100000);
          RayRenderVRML1(ray, ray_width, ray_height, &vla,
                         I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov, angle, I->m_view.m_pos[2]);
          *charVLA_ptr = vla;
        }
        break;
      case cSceneRay_MODE_IDTF:
        {
          *headerVLA_ptr = VLACalloc(char, 10000);
          *charVLA_ptr = VLACalloc(char, 10000);
          RayRenderIDTF(ray, headerVLA_ptr, charVLA_ptr);
        }
        break;
      case 8:                   /* mode 8 is COLLADA (.dae) */
        {
          *charVLA_ptr = VLACalloc(char, 100000);
          RayRenderCOLLADA(ray, ray_width, ray_height, charVLA_ptr,
                            I->m_view.m_clipSafe.m_front, I->m_view.m_clipSafe.m_back, fov);
        }
        break;

      }
      RayFree(ray);
    }
    if(I->grid.active)
      GridSetRayViewport(&I->grid, -1, &ray_x, &ray_y, &ray_width, &ray_height);

    if((mode == 0) && I->Image && !I->Image->empty()) {
      SceneApplyImageGamma(G, I->Image->pixels(), I->Image->getWidth(),
                           I->Image->getHeight());
    }

    stereo_hand--;
    if((I->StereoMode == 0) || (stereo_hand <= 0))
      break;
    else {
      stereo_image = I->Image;
    }
  }

  if(stereo_image) {
    if(I->Image) {
      switch (I->StereoMode) {
      case cStereo_quadbuffer:
      case cStereo_geowall:
      case cStereo_openvr:
        /* merge the two images into one */
        I->Image->merge(*stereo_image);
        break;
      case cStereo_crosseye:
      case cStereo_walleye:
        {
          /* merge the two images into one */
          auto merged_image =
              pymol::Image(I->Image->getWidth() * 2, I->Image->getHeight());

          unsigned int *q = merged_image.pixels();
          unsigned int *l;
          unsigned int *r;
          int height, width;
          int a, b;

          if(I->StereoMode == 2) {
            l = (unsigned int *) stereo_image->bits();
            r = (unsigned int *) I->Image->bits();
          } else {
            r = (unsigned int *) stereo_image->bits();
            l = (unsigned int *) I->Image->bits();
          }
          height = I->Image->getHeight();
          width = I->Image->getWidth();

          for(a = 0; a < height; a++) {
            for(b = 0; b < width; b++)
              *(q++) = *(l++);
            for(b = 0; b < width; b++)
              *(q++) = *(r++);
          }
          *I->Image = std::move(merged_image);
        }
        break;
      case cStereo_anaglyph:
        {
          int big_endian;
          {
            unsigned int test;
            unsigned char *testPtr;
            test = 0xFF000000;
            testPtr = (unsigned char *) &test;
            big_endian = (*testPtr) & 0x01;
          }
          {
            extern float anaglyphR_constants[6][9];
            extern float anaglyphL_constants[6][9];
            unsigned int *l = stereo_image->pixels();
            unsigned int *r = I->Image->pixels();
	    int anaglyph_mode = SettingGetGlobal_i(G, cSetting_anaglyph_mode);
	    /* anaglyph scalars */
	    float * a_r = anaglyphR_constants[anaglyph_mode];
	    float * a_l = anaglyphL_constants[anaglyph_mode];

            int height, width;
            int a, b;
	    float _r[3] = {0.F,0.F,0.F}, _l[3] = {0.F,0.F,0.F}, _b[3] = {0.F,0.F,0.F};
            height = I->Image->getHeight();
            width = I->Image->getWidth();
            
            for(a = 0; a < height; a++) {
              for(b = 0; b < width; b++) {
                if(big_endian) {
                  /* original : RGBA
		   *r = (*l & 0x00FFFFFF) | (*r & 0xFF000000);
		   */
		  /* UNTESTED */
		  _l[0] = (float)((*r & 0xFF000000));
		  _l[1] = (float)((*r & 0x00FF0000) >> 16);
		  _l[2] = (float)((*r & 0x0000FF00) >> 8);
		  _r[0] = (float)((*l & 0xFF000000));
		  _r[1] = (float)((*l & 0x00FF0000) >> 16);
		  _r[2] = (float)((*l & 0x0000FF00) >> 8);
		  _b[0] = (a_l[0] * _l[0] + a_l[3] * _l[1] + a_l[6] * _l[2]); // R
		  _b[1] = (a_l[1] * _l[0] + a_l[4] * _l[1] + a_l[7] * _l[2]); // G
		  _b[2] = (a_l[2] * _l[0] + a_l[5] * _l[1] + a_l[8] * _l[2]); // B
		  *l = (unsigned int) (0x000000FF & *l) |
		       (unsigned int) (1.0 * _b[0]) |
		      ((unsigned int) (1.0 * _b[1]))<<8 |
		      ((unsigned int) (1.0 * _b[2]))<<16;

		  _b[0] = (a_r[0] * _r[0] + a_r[3] * _r[1] + a_r[6] * _r[2]); // R
		  _b[1] = (a_r[1] * _r[0] + a_r[4] * _r[1] + a_r[7] * _r[2]); // G
		  _b[2] = (a_r[2] * _r[0] + a_r[5] * _r[1] + a_r[8] * _r[2]); // B

		  *r = (unsigned int) (0x000000FF & *r) |
		       (unsigned int) (1.0 * _b[0]) |
   		      ((unsigned int) (1.0 * _b[1]))<<8 |
		      ((unsigned int) (1.0 * _b[2]))<<16;

		  *r = (*l | *r);
                } else {
                  /* original : AGBR
		   *r = (*l & 0xFFFFFF00) | (*r & 0x000000FF);
		   */
		  
		  /* Right and Left as unsigned ints */
		  /* CORRECT */
		  _l[0] = (float)((*r & 0x000000FF));
		  _l[1] = (float)((*r & 0x0000FF00) >> 8);
		  _l[2] = (float)((*r & 0x00FF0000) >> 16);
		  _r[0] = (float)((*l & 0x000000FF));
		  _r[1] = (float)((*l & 0x0000FF00) >> 8);
		  _r[2] = (float)((*l & 0x00FF0000) >> 16);

		  _b[0] = (a_l[0] * _l[0] + a_l[3] * _l[1] + a_l[6] * _l[2]); // R
		  _b[1] = (a_l[1] * _l[0] + a_l[4] * _l[1] + a_l[7] * _l[2]); // G
		  _b[2] = (a_l[2] * _l[0] + a_l[5] * _l[1] + a_l[8] * _l[2]); // B

		  *l = (unsigned int) (0xFF000000 & *l) |
		       (unsigned int) (1.0 * _b[0]) |
		      ((unsigned int) (1.0 * _b[1]))<<8 |
		      ((unsigned int) (1.0 * _b[2]))<<16;

		  _b[0] = (a_r[0] * _r[0] + a_r[3] * _r[1] + a_r[6] * _r[2]); // R
		  _b[1] = (a_r[1] * _r[0] + a_r[4] * _r[1] + a_r[7] * _r[2]); // G
		  _b[2] = (a_r[2] * _r[0] + a_r[5] * _r[1] + a_r[8] * _r[2]); // B

		  *r = (unsigned int) (0xFF000000 & *r) |
		       (unsigned int) (1.0 * _b[0]) |
   		      ((unsigned int) (1.0 * _b[1]))<<8 |
		      ((unsigned int) (1.0 * _b[2]))<<16;

		  *r = (*l | *r);
                }
                l++;
                r++;
              }
            }
          }
        }
        break;
      case cStereo_stencil_by_row:
      case cStereo_stencil_by_column:
      case cStereo_stencil_checkerboard:
        {
          /* merge the two images into one */

          int parity = 0;

          if(I->StereoMode == cStereo_stencil_by_row) {
            parity = I->StencilParity;
            if(I->rect.bottom & 0x1)
              parity = 1 - parity;
          }

          unsigned int* l = stereo_image->pixels();
          unsigned int* r = I->Image->pixels();

          int height = I->Image->getHeight();
          int width = I->Image->getWidth();

          auto merged_image = pymol::Image(width, height);
          unsigned int *q = merged_image.pixels();

          for (int a = 0; a < height; ++a) {
            for (int b = 0; b < width; ++b) {
              switch (I->StereoMode) {
              case cStereo_stencil_by_row:
                if((a + parity) & 0x1) {
                  *(q++) = *(l++);
                  r++;
                } else {
                  *(q++) = *(r++);
                  l++;
                }
                break;
              case cStereo_stencil_by_column:
                if(b & 0x1) {
                  *(q++) = *(l++);
                  r++;
                } else {
                  *(q++) = *(r++);
                  l++;
                }
                break;
              case cStereo_stencil_checkerboard:
                if((a + b) & 0x1) {
                  *(q++) = *(l++);
                  r++;
                } else {
                  *(q++) = *(r++);
                  l++;
                }
                break;
              }
            }
          }
          *I->Image = std::move(merged_image);
        }   
        break;
      }
    }
  }
  timing = UtilGetSeconds(G) - timing;
  if(mode != 2) {               /* don't show timings for tests */
    accumTiming += timing;

    if(show_timing && !quiet) {
      if(!G->Interrupt) {
        PRINTFB(G, FB_Ray, FB_Details)
          " Ray: render time: %4.2f sec. = %3.1f frames/hour (%4.2f sec. accum.).\n",
          timing, 3600 / timing, accumTiming ENDFB(G);
      } else {
        PRINTFB(G, FB_Ray, FB_Details)
          " Ray: render aborted.\n" ENDFB(G);
      }
    }
  }

  if(mode != 3) {
    OrthoDirty(G);
  }

  /* EXPERIMENTAL VOLUME CODE */
  if (rayVolume) {
    SceneUpdate(G, true);
  }
  OrthoBusyFast(G, 20, 20);
  PyMOL_SetBusy(G->PyMOL, false);

  return true;
#endif
}

static int SceneDeferredRay(DeferredRay * dr)
{
  PyMOLGlobals *G = dr->m_G;
  SceneRay(G, dr->ray_width, dr->ray_height, dr->mode,
           NULL, NULL, dr->angle, dr->shift, dr->quiet,
           NULL, dr->show_timing, dr->antialias);
  if((dr->mode == 0) && G->HaveGUI && SettingGetGlobal_b(G, cSetting_auto_copy_images)) {
#ifdef _PYMOL_IP_EXTRAS
    PParse(G, "cmd._copy_image(quiet=0)");
#else
#ifdef PYMOL_EVAL
    PRINTFB(G, FB_Scene, FB_Warnings)
      " Warning: Clipboard image transfers disabled in Evaluation Builds.\n" ENDFB(G);
#endif
#endif
  }
  return 1;
}

int SceneDeferRay(PyMOLGlobals * G,
                  int ray_width,
                  int ray_height,
                  int mode,
                  float angle, float shift, int quiet, int show_timing, int antialias)
{
  auto dr = pymol::make_unique<DeferredRay>(G);
  if(dr) {
    dr->ray_width = ray_width;
    dr->ray_height = ray_height;
    dr->mode = mode;
    dr->angle = angle;
    dr->shift = shift;
    dr->quiet = quiet;
    dr->show_timing = show_timing;
    dr->antialias = antialias;
    dr->fn = (DeferredFn *) SceneDeferredRay;
  }
  OrthoDefer(G, std::move(dr));
  return 1;
}

void SceneRenderRayVolume(PyMOLGlobals * G, CScene *I){
#ifndef PURE_OPENGL_ES_2
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, I->Width, 0, I->Height, -100, 100);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
#endif
  
#ifndef PURE_OPENGL_ES_2
  glRasterPos3f(0, 0, -1);
#endif
  glDepthMask(GL_FALSE);
#ifndef PURE_OPENGL_ES_2
  if (PIsGlutThread() && rayVolumeImage) {
    if (rayWidth == I->Width && rayHeight == I->Height){
      glDrawPixels(rayVolumeImage->getWidth(), rayVolumeImage->getHeight(),
          GL_RGBA, GL_UNSIGNED_BYTE, rayVolumeImage->bits());
    } else {
      SceneDrawImageOverlay(G, 1, NULL);
    }
  }
#endif
  glDepthMask(GL_TRUE);
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
  glDepthFunc(GL_ALWAYS);
#ifndef PURE_OPENGL_ES_2
  if (PIsGlutThread() && rayWidth == I->Width && rayHeight == I->Height)
    glDrawPixels(I->Width, I->Height, GL_DEPTH_COMPONENT, GL_FLOAT, rayDepthPixels); 
#endif
  glDepthFunc(GL_LESS);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  
#ifdef PURE_OPENGL_ES_2
  /* TODO */
#else
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
#endif
}
