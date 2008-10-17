/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* TJO: TJ O'Donnell, under contract for DeLano Scientific LLC
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_predef.h"

/* backwards compatibility */

#ifdef _HAVE_LIBPNG
#ifndef _PYMOL_LIBPNG
#define _PYMOL_LIBPNG
#endif
#endif

/* do we have libpng? */

#ifdef _PYMOL_LIBPNG
#include<png.h>

 /* The png_jmpbuf() macro, used in error handling, became available in
  * libpng version 1.0.6.  If you want to be able to run your code with older
  * versions of libpng, you must define the macro yourself (but only if it
  * is not already defined by libpng!).
  */

#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif

#endif

#include"os_std.h"

#include"Base.h"
#include "MyPNG.h"
#include"MemoryDebug.h"
#include "Setting.h"


int MyPNGWrite(PyMOLGlobals *G,char *file_name,unsigned char *p,
               unsigned int width,unsigned int height,float dpi)
{
#ifdef _PYMOL_LIBPNG
  int ok=true;
  FILE *fp = NULL;
  png_structp png_ptr;
  png_infop info_ptr;
  int bit_depth = 8;
  int bytes_per_pixel = 4;
  png_uint_32 k;
  png_byte *image = (png_byte*)p;
  png_bytep *row_pointers;
  int fd = 0;

  row_pointers=Alloc(png_bytep,height);
  
  /* open the file, allowing use of an encoded file descriptor, with
     approach adapted from TJO: chr(1) followed by ascii-format integer */
  if(file_name[0] == 1) {
    if( sscanf(file_name+1, "%d", &fd) == 1) {
      fp = fdopen(fd, "wb");
    }
  } else {
    fp = fopen(file_name, "wb");
  }
  if (fp == NULL) {
    ok=false;
    goto cleanup;
  } else if(feof(fp)) {
    ok=false;
    goto cleanup;
  }
   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);

   if (png_ptr == NULL) {
     ok=false;
     goto cleanup;
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
     png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
     ok=false;
     goto cleanup;
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr))) {
     /* If we get here, we had a problem reading the file */
     png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
     ok=false;
     goto cleanup;
   }

   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, PNG_COLOR_TYPE_RGB_ALPHA,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   if(dpi>0.0F) { /* only set resolution if dpi is positive */
     int dots_per_meter = (int)(dpi*39.3700787);
     png_set_pHYs(png_ptr, info_ptr, dots_per_meter, dots_per_meter, PNG_RESOLUTION_METER);
   }


   png_set_gamma(png_ptr, SettingGet(G,cSetting_png_screen_gamma), 
                 SettingGet(G,cSetting_png_file_gamma));
   
   /* stamp the image as being created by PyMOL we could consider
    * supporting optional annotations as well: PDB codes, canonical
    * smiles, INCHIs, and other common identifiers */

   {
     png_text text;
     text.compression = PNG_TEXT_COMPRESSION_NONE;
     text.key = (png_charp) "Software";
     text.text = (png_charp) "PyMOL";
     text.text_length = 5;
     png_set_text(png_ptr, info_ptr, &text, 1);
   }
   {
     png_text text;
     text.compression = PNG_TEXT_COMPRESSION_NONE;
     text.key = (png_charp) "URL";
     text.text = (png_charp) "http://www.pymol.org";
     text.text_length = 5;
     png_set_text(png_ptr, info_ptr, &text, 1);
   }

   /* Write the file header information.  REQUIRED */
   png_write_info(png_ptr, info_ptr);

   /* The easiest way to write the image (you may have a different memory
    * layout, however, so choose what fits your needs best).  You need to
    * use the first method if you aren't handling interlacing yourself.
    */
   for (k = 0; k < height; k++)
     row_pointers[(height-k)-1] = image + k*width*bytes_per_pixel;

   png_write_image(png_ptr, row_pointers);

   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

 cleanup:
   if(fp) {
     /* close the file */
     fclose(fp);
   }
   
   mfree(row_pointers);
   /* that's it */

   return ok;
#else

   return 0;
#endif
}

int MyPNGRead(char *file_name,unsigned char **p_ptr,unsigned int *width_ptr,unsigned int *height_ptr)
{

#ifdef _PYMOL_LIBPNG

  FILE *png_file=NULL;
  png_struct    *png_ptr = NULL;
  png_info	*info_ptr = NULL;
  png_byte      buf[8];
  png_byte      *png_pixels = NULL;
  png_byte      **row_pointers = NULL;
  png_byte      *pix_ptr = NULL;
  png_uint_32   row_bytes=0;

  png_uint_32   width;
  png_uint_32   height;
  int           bit_depth;
  int           color_type;
  int           row, col;
  int           ret;
  int           i;
  int ok=true;
  unsigned char *p=NULL;
  double        file_gamma;

  if(!file_name)
    return 0;

   png_file = fopen(file_name, "rb");
   if (png_file == NULL)
     return 0;

   /* read and check signature in PNG file */
   ret = fread (buf, 1, 8, png_file);
   if (ret != 8)
     ok=false;
   
   if(ok) {
     ret = png_check_sig (buf, 8);
     if (!ret)
       ok=false;
   }
   /* create png and info structures */
   if(ok) {
     png_ptr = png_create_read_struct (PNG_LIBPNG_VER_STRING,
                                       NULL, NULL, NULL);
     if (!png_ptr)
       ok=false;
   }
   
   if(ok) {
     info_ptr = png_create_info_struct (png_ptr);
     if (!info_ptr)
       ok=false;
   }
   
   if (setjmp (png_jmpbuf(png_ptr)))
     ok = false;
   
   if(ok) {
     /* set up the input control for C streams */
     png_init_io (png_ptr, png_file);
     png_set_sig_bytes (png_ptr, 8);  /* we already read the 8 signature bytes */
     
     /* read the file information */
     png_read_info (png_ptr, info_ptr);
     
     /* get size and bit-depth of the PNG-image */
     png_get_IHDR (png_ptr, info_ptr,
                   &width, &height, &bit_depth, &color_type,
                   NULL, NULL, NULL);
     
     /* set-up the transformations */

     if(color_type!=PNG_COLOR_TYPE_RGB_ALPHA) {
       png_set_expand(png_ptr);
       png_set_filler(png_ptr,0xFF,PNG_FILLER_AFTER);
     }
     
     if (color_type == PNG_COLOR_TYPE_GRAY ||
         color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
       png_set_gray_to_rgb (png_ptr);
     /* only if file has a file gamma, we do a correction */
     if (png_get_gAMA (png_ptr, info_ptr, &file_gamma))
       png_set_gamma (png_ptr, (double) 2.2, file_gamma);
     
     /* all transformations have been registered; now update info_ptr data,
      * get rowbytes and channels, and allocate image memory */
     
     png_read_update_info (png_ptr, info_ptr);
     
     /* get the new color-type and bit-depth (after expansion/stripping) */
     png_get_IHDR (png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
                   NULL, NULL, NULL);
     
     /* row_bytes is the width x number of channels x (bit-depth / 8) */
     row_bytes = png_get_rowbytes (png_ptr, info_ptr);
     if ((png_pixels = (png_byte *) malloc (row_bytes * height * sizeof (png_byte))) == NULL) {
       ok=false;
     }
   }
   
   if(ok) {
     
     if ((row_pointers = (png_byte **) malloc (height * sizeof (png_bytep))) == NULL)
       {
         png_destroy_read_struct (&png_ptr, &info_ptr, NULL);
         free (png_pixels);
         png_pixels = NULL;
         ok=false;
       }
   }
   
   if(ok) {
     
     /* set the individual row_pointers to point at the correct offsets */
     for (i = 0; i < ((signed)height); i++)
       row_pointers[i] = png_pixels + i * row_bytes;
     
     /* now we can go ahead and just read the whole image */
     png_read_image (png_ptr, row_pointers);
     
     /* read rest of file, and get additional chunks in info_ptr - REQUIRED */
     png_read_end (png_ptr, info_ptr);
   }
   
   if(ok) {
     /* now reformat image into PyMOL format */
     
     p=(unsigned char*)mmalloc(4*width*height);
     if(!p)
       ok=false;
   }
   if(ok) {
     *(p_ptr)=p;
     *(width_ptr)=width;
     *(height_ptr)=height;
     
     for (row = 0; row < (signed)height; row++)
       {
         pix_ptr=row_pointers[(height-1)-row];
         for (col = 0; col < (signed)width; col++)
           {
             *p++=*pix_ptr++;
             *p++=*pix_ptr++;
             *p++=*pix_ptr++;
             *p++=*pix_ptr++;
           }
       }
     
   } 
   
   if (row_pointers != (unsigned char**) NULL)
     free (row_pointers);
   if (png_pixels != (unsigned char*) NULL)
     free (png_pixels);
   
   if(png_ptr) {
     png_destroy_read_struct (&png_ptr, &info_ptr, (png_infopp) NULL);
   }
   if(png_file)
     fclose(png_file);

   return(ok);
#else
   return (false);
#endif
  
} /* end of source */









