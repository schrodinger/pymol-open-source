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
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_std.h"

#include "MyPNG.h"
#include"MemoryDebug.h"

#ifdef _HAVE_LIBPNG
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


void MyPNGWrite(char *file_name,unsigned char *p,unsigned int width,unsigned int height)
{
#ifdef _HAVE_LIBPNG

   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
	int bit_depth = 8;
	int bytes_per_pixel = 4;
   png_uint_32 k;
   png_byte *image = (png_byte*)p;
   png_bytep *row_pointers;

	row_pointers=Alloc(png_bytep,height);

   /* open the file */
   fp = fopen(file_name, "wb");
   if (fp == NULL)
      return;

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);

   if (png_ptr == NULL)
   {
      fclose(fp);
      return;
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
      return;
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
      return;
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

   png_set_gamma(png_ptr, 2.4, 1.0);

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

   /* close the file */
   fclose(fp);

	mfree(row_pointers);
   /* that's it */

#endif

   return;
}


