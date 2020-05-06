/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* TJO: TJ O'Donnell, under contract for Schrodinger, LLC
-* 
-*
Z* -------------------------------------------------------------------
*/
#include"os_predef.h"
#include"MemoryDebug.h"
#include"pymol/memory.h"


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
#include "Err.h"
#include "File.h"

/*
 * base64 decoding
 * http://stackoverflow.com/questions/342409/how-do-i-base64-encode-decode-in-c
 */

static unsigned char base64_decoding_table[] = {
  62,   0,   0,   0,  63,  52,  53,  54,
  55,  56,  57,  58,  59,  60,  61,   0,
   0,   0,   0,   0,   0,   0,   0,   1,
   2,   3,   4,   5,   6,   7,   8,   9,
  10,  11,  12,  13,  14,  15,  16,  17,
  18,  19,  20,  21,  22,  23,  24,  25,
   0,   0,   0,   0,   0,   0,  26,  27,
  28,  29,  30,  31,  32,  33,  34,  35,
  36,  37,  38,  39,  40,  41,  42,  43,
  44,  45,  46,  47,  48,  49,  50,  51};

static unsigned char *base64_decode(const char *data,
				    size_t input_length=0) {
  unsigned char *decoded_data = NULL;
  unsigned int triple;
  unsigned int i = 0, j = 0, k;
  char c;

  if (input_length < 1)
    input_length = strlen(data);

  ok_assert(1, decoded_data = pymol::malloc<unsigned char>(input_length / 4 * 3));

  while (i < input_length) {
    triple = 0;

    for (k = 4; k && i < input_length;) {
      c = data[i++];
      if (c < 43 || c > 122)
        continue;
      triple += base64_decoding_table[c - 43] << (--k) * 6;
    }

    ok_assert(1, k == 0);

    for (k = 3; k;)
      decoded_data[j++] = (triple >> (--k) * 8) & 0xFF;
  }

  return decoded_data;

ok_except1:
  mfree(decoded_data);
  return NULL;
}

typedef struct {
  unsigned char *c; // moving pointer
  unsigned char *h; // for freeing
} uchar2p;

#ifdef _PYMOL_LIBPNG
/**
 * Use with png_set_read_fn, for reading a PNG file from memory.
 * This simply copies from the memory pointer to outBytes while
 * incrementing the memory pointer.
 */
void read_data_from_buffer(png_structp png_ptr,
                           png_bytep outBytes,
                           png_size_t byteCountToRead) {
  uchar2p *io_ptr = (uchar2p*) png_get_io_ptr(png_ptr);

  if(!io_ptr)
    return;

  while(byteCountToRead--) {
    *(outBytes++) = *(io_ptr->c++);
  }
}

/**
 * Use with png_set_write_fn, for writing a PNG file to memory.
 */
static void write_data_to_buffer(png_structp png_ptr,
                                 png_bytep data,
                                 png_size_t length) {
  auto io_ptr = reinterpret_cast<png_outbuf_t*>(png_get_io_ptr(png_ptr));
  io_ptr->insert(io_ptr->end(), data, data + length);
}

/**
 * Use with png_set_read_fn instead of png_init_io, allows mixing of
 * Visual Studio versions
 */
static void read_data_from_file(
    png_structp png_ptr,
    png_bytep buffer,
    png_size_t count) {
  auto fp = static_cast<FILE*>(png_get_io_ptr(png_ptr));
  fread(buffer, 1, count, fp);
}

/**
 * Use with png_set_write_fn instead of png_init_io, allows mixing of
 * Visual Studio versions
 */
static void write_data_to_file(
    png_structp png_ptr,
    png_bytep buffer,
    png_size_t count) {
  auto fp = static_cast<FILE*>(png_get_io_ptr(png_ptr));
  fwrite(buffer, 1, count, fp);
}
#endif

int MyPNGWrite(const char* file_name, const pymol::Image& img, const float dpi,
    const int format, const int quiet, const float screen_gamma,
    const float file_gamma, png_outbuf_t* io_ptr)
{
  const unsigned char* data_ptr = img.bits();
  int width = img.getWidth();
  int height = img.getHeight();
  switch (format) {
  case cMyPNG_FormatPNG:
    {
#ifdef _PYMOL_LIBPNG
      int ok = true;
      FILE *fp = NULL;
      png_structp png_ptr;
      png_infop info_ptr;
      int bit_depth = 8;
      int bytes_per_pixel = 4;
      png_uint_32 k;
      png_byte *image = (png_byte *) data_ptr;
      png_bytep *row_pointers;
      int fd = 0;

      row_pointers = pymol::malloc<png_bytep>(height);

      /* open the file, allowing use of an encoded file descriptor, with
         approach adapted from TJO: chr(1) followed by ascii-format integer */
      if (!io_ptr) {
        if(file_name[0] == 1) {
          if(sscanf(file_name + 1, "%d", &fd) == 1) {
            fp = fdopen(fd, "wb");
          }
        } else {
          fp = pymol_fopen(file_name, "wb");
        }
        if(fp == NULL) {
          ok = false;
          goto cleanup;
        } else if(feof(fp)) {
          ok = false;
          goto cleanup;
        }
      }
      /* Create and initialize the png_struct with the desired error handler
       * functions.  If you want to use the default stderr and longjump method,
       * you can supply NULL for the last three parameters.  We also check that
       * the library version is compatible with the one used at compile time,
       * in case we are using dynamically linked libraries.  REQUIRED.
       */
      png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

      if(png_ptr == NULL) {
        ok = false;
        goto cleanup;
      }

      /* Allocate/initialize the image information data.  REQUIRED */
      info_ptr = png_create_info_struct(png_ptr);
      if(info_ptr == NULL) {
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        ok = false;
        goto cleanup;
      }

      /* Set error handling.  REQUIRED if you aren't supplying your own
       * error handling functions in the png_create_write_struct() call.
       */
      if(setjmp(png_jmpbuf(png_ptr))) {
        /* If we get here, we had a problem reading the file */
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
        ok = false;
        goto cleanup;
      }

      if (io_ptr) {
        png_set_write_fn(png_ptr, (void*) io_ptr, write_data_to_buffer, NULL);
      } else {
        /* set up the output control if you are using standard C streams */
        png_set_write_fn(png_ptr, (void*) fp, write_data_to_file, NULL);
      }

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

      if(dpi > 0.0F) {          /* only set resolution if dpi is positive */
        int dots_per_meter = (int) (dpi * 39.3700787);
        png_set_pHYs(png_ptr, info_ptr, dots_per_meter, dots_per_meter,
                     PNG_RESOLUTION_METER);
      }

      png_set_gamma(png_ptr, screen_gamma, file_gamma);

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
      for(k = 0; k < height; k++)
        row_pointers[(height - k) - 1] = image + k * width * bytes_per_pixel;

      png_write_image(png_ptr, row_pointers);

      /* It is REQUIRED to call this to finish writing the rest of the file */
      png_write_end(png_ptr, info_ptr);

      /* clean up after the write, and free any memory allocated */
      png_destroy_write_struct(&png_ptr, &info_ptr);

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
    break;
  case cMyPNG_FormatPPM:
    {
      FILE *fil = pymol_fopen(file_name, "wb");
      unsigned char *buffer = pymol::malloc<unsigned char>(3 * width * height);

      if(fil && buffer) {
        fprintf(fil, "P6\n");
        fprintf(fil, "%d %d\n", width, height);
        fprintf(fil, "255\n");
        {
          unsigned int a, b;
          unsigned char *q = buffer;
          const unsigned char *p;
          p = data_ptr + width * 4 * (height - 1);
          for(b = 0; b < height; b++) {
            for(a = 0; a < width; a++) {
              *(q++) = *(p++);  /* platform-specific ordering? */
              *(q++) = *(p++);
              *(q++) = *(p++);
              p++;
            }
            p -= width * 8;
          }
          fwrite(buffer, width, height * 3, fil);
        }
      }
      if(fil) {
        fclose(fil);
      }
      FreeP(buffer);
    }
    return 1;
    break;
  }
  return 0;
}

std::unique_ptr<pymol::Image> MyPNGRead(const char *file_name)
{
  std::unique_ptr<pymol::Image> img;
#ifdef _PYMOL_LIBPNG

  FILE *png_file = NULL;
  png_struct *png_ptr = NULL;
  png_info *info_ptr = NULL;
  png_byte buf[8];
  png_byte *png_pixels = NULL;
  png_byte **row_pointers = NULL;
  png_byte *pix_ptr = NULL;
  png_uint_32 row_bytes = 0;

  png_uint_32 width;
  png_uint_32 height;
  int bit_depth;
  int color_type;
  int row, col;
  int ret;
  int i;
  int ok = true;
  double file_gamma;
  uchar2p data = {NULL, NULL};

  if(!file_name)
    return nullptr;

  if(!strncmp(file_name, "data:image/png;base64,", 22)) {
    const char *base64str = file_name + 22;
    data.h = data.c = base64_decode(base64str);
    memcpy(buf, data.c, 8);
    data.c += 8;

  } else {
    png_file = pymol_fopen(file_name, "rb");
    if(png_file == NULL)
      return nullptr;

    /* read and check signature in PNG file */
    ret = fread(buf, 1, 8, png_file);
    if(ret != 8)
      ok = false;
  }

  if(ok) {
    ret = png_check_sig(buf, 8);
    if(!ret)
      ok = false;
  }
  /* create png and info structures */
  if(ok) {
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_ptr)
      ok = false;
  }

  if(ok) {
    info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr)
      ok = false;
  }

  if(ok && setjmp(png_jmpbuf(png_ptr)))
    ok = false;

  if(ok) {
    /* set up the input control for C streams */

    if(data.h) {
      png_set_read_fn(png_ptr, (void*) &data, read_data_from_buffer);
    } else {
      png_set_read_fn(png_ptr, (void*) png_file, read_data_from_file);
    }

    png_set_sig_bytes(png_ptr, 8);      /* we already read the 8 signature bytes */

    /* read the file information */
    png_read_info(png_ptr, info_ptr);

    /* get size and bit-depth of the PNG-image */
    png_get_IHDR(png_ptr, info_ptr,
                 &width, &height, &bit_depth, &color_type, NULL, NULL, NULL);

    /* set-up the transformations */

    if(color_type != PNG_COLOR_TYPE_RGB_ALPHA) {
      png_set_expand(png_ptr);
      png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);
    }

    if(color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
      png_set_gray_to_rgb(png_ptr);
    /* only if file has a file gamma, we do a correction */
    if(png_get_gAMA(png_ptr, info_ptr, &file_gamma))
      png_set_gamma(png_ptr, (double) 2.2, file_gamma);

    /* all transformations have been registered; now update info_ptr data,
     * get rowbytes and channels, and allocate image memory */

    png_read_update_info(png_ptr, info_ptr);

    /* get the new color-type and bit-depth (after expansion/stripping) */
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
                 NULL, NULL, NULL);

    /* row_bytes is the width x number of channels x (bit-depth / 8) */
    row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    if((png_pixels = (png_byte *) malloc(row_bytes * height * sizeof(png_byte))) == NULL) {
      ok = false;
    }
  }

  if(ok) {

    if((row_pointers = (png_byte **) malloc(height * sizeof(png_bytep))) == NULL) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      free(png_pixels);
      png_pixels = NULL;
      ok = false;
    }
  }

  if(ok) {

    /* set the individual row_pointers to point at the correct offsets */
    for(i = 0; i < ((signed) height); i++)
      row_pointers[i] = png_pixels + i * row_bytes;

    /* now we can go ahead and just read the whole image */
    png_read_image(png_ptr, row_pointers);

    /* read rest of file, and get additional chunks in info_ptr - REQUIRED */
    png_read_end(png_ptr, info_ptr);
  }

  img = pymol::make_unique<pymol::Image>(width, height);
  if(ok) {
    auto p = img->bits();
    for(row = 0; row < (signed) height; row++) {
      pix_ptr = row_pointers[(height - 1) - row];
      for(col = 0; col < (signed) width; col++) {
        *p++ = *pix_ptr++;
        *p++ = *pix_ptr++;
        *p++ = *pix_ptr++;
        *p++ = *pix_ptr++;
      }
    }

  }

  if(row_pointers != (unsigned char **) NULL)
    free(row_pointers);
  if(png_pixels != (unsigned char *) NULL)
    free(png_pixels);

  if(png_ptr) {
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);
  }
  if(png_file)
    fclose(png_file);
  if(data.h)
    mfree(data.h);

  return img;
#else
  return nullptr;
#endif

}                               /* end of source */
