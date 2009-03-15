/* MACHINE GENERATED FILE, DO NOT EDIT! */

#define VMDPLUGIN molfile_dtrplugin
#define STATIC_PLUGIN 1

//
// Version info for VMD plugin tree:
//   $Id: dtrplugin.cxx,v 1.15 2009/03/02 15:59:47 johns Exp $
//
// Version info for last sync with D. E. Shaw Research:
//  //depot/desrad/main/sw/libs/vmd_plugins,DIST/dtrplugin.cxx#3
//

/*
Copyright 2009, D. E. Shaw Research, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research, LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#if defined(_MSC_VER)
#ifndef DESRES_WIN32
#define DESRES_WIN32
#endif
#endif

#include <stdio.h>
#ifdef DESRES_WIN32
#include <io.h>
#include <direct.h>
#include <fcntl.h>
#include <windows.h>
typedef int int32_t;
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
#if 1
typedef unsigned __int64 uint64_t;    // This also works with MVSC6
#else
typedef unsigned long long uint64_t;
#endif
typedef unsigned short uint16_t;
typedef unsigned int ssize_t;
typedef int mode_t;
#define mkdir(a,b) _mkdir(a)
#define rmdir(a)   _rmdir(a)
#define ftello(a)  ftell(a)
#else
#define O_BINARY 0
#include <inttypes.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#endif
#include <sstream>
#include <ios>
#include <iomanip>
#include <math.h>
#include <errno.h>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>
#include <ios>
#include <set>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "vmddir.h"
#include "endianswap.h"

#ifndef DESRES_WIN32
static const char s_sep = '/';
#include <sys/mman.h>

#include <netinet/in.h> /* for htonl */
#if defined(_AIX)
#include <fcntl.h>
#else
#include <sys/fcntl.h>
#endif

#if defined(__sun)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#else
/// windows version

#define M_PI (3.1415926535897932385)
#define M_PI_2 (1.5707963267948966192)

#ifndef S_ISDIR
#define S_ISDIR(m) (((m) & _S_IFDIR) == _S_IFDIR)
#endif

static const char s_sep = '\\';
#define PROT_READ (0)
#define PROT_WRITE (0)
#define MAP_PRIVATE (0)
#define MAP_FAILED ((void*)-1)
static void * mmap( void *start, size_t length, int prot, int flags, 
                    int fd, off_t offset) {
  ::_lseek(fd,offset,SEEK_SET);
  start = ::malloc(length);
  // TODO: interrupting read problem
  ::read(fd,reinterpret_cast<char*>(start),length);
  return start;
}

static int munmap(void *start, size_t length) {
  ::free(start);
  return 0;
}
#endif

#if defined(__sun)
#define MMAPTYPE(a) ((caddr_t) a)
#else
#define MMAPTYPE(a) (a)
#endif

namespace {

  /// definitions of binary representation of frameset files

  typedef struct {
    uint32_t magic;           //!< Magic number
    uint32_t version;         //!< Version of creator
    uint32_t framesize_lo;    //!< bytes in frame (low)
    uint32_t framesize_hi;    //!< bytes in frame (high)
  } required_header_t;

  //! Header structure within file.
  typedef struct {
    required_header_t required; //!< 4 word mini-header

    uint32_t size_header_block; //!< Size of this header
    uint32_t unused0;         //!< not used in current implementation
    uint32_t irosetta;        //!< 32-bit integer rosetta value.
    float    frosetta;        //!< 32-bit float rosetta

    uint32_t drosetta_lo;     //!< 64-bit float rosetta (low)
    uint32_t drosetta_hi;     //!< 64-bit float rosetta (high)
    uint32_t lrosetta_lo;     //!< 64-bit integer rosetta (low)
    uint32_t lrosetta_hi;     //!< 64-bit integer rosetta (high)

    uint32_t endianism;       //!< Endianism of writer machine.
    uint32_t nlabels;         //!< Number of labeled fields.
    uint32_t size_meta_block; //!< Number of bytes of meta information (padded)
    uint32_t size_typename_block; //!< Number of bytes of typenames (padded)

    uint32_t size_label_block; //!< Number of bytes to store label strings (padded)
    uint32_t size_scalar_block; //!< Number of bytes of scalar storage (padded)
    uint32_t size_field_block_lo; //!< Number of bytes of field storage (padded)
    uint32_t size_field_block_hi; //!< Number of bytes of field storage (padded)

    uint32_t size_crc_block;  //!< Size of the trailing CRC field (unused!)
    uint32_t size_padding_block; //!< Number of ignored bytes to pad to pagesize boundary.
    uint32_t unused1;         //!< Not used in current implementation.
    uint32_t unused2;         //!< Not used in current implementation.

  } header_t;

  typedef struct {
    uint32_t type;            //!< \brief Typecode for this type.
    uint32_t elementsize;     //!< \brief Number of bytes in an element
    uint32_t count_lo;        //!< \brief Number of elements (low)
    uint32_t count_hi;        //!< \brief Number of elements (high)
  } metadisk_t;

  typedef struct key_prologue {
    uint32_t magic;           /* Magic number for frames */
    uint32_t frames_per_file; /* Number of frames in each file */
    uint32_t key_record_size; /* The size of each key record */
  } key_prologue_t;

  typedef struct key_record {
    uint32_t time_lo;       /* Time associated with this file (low bytes). */
    uint32_t time_hi;       /* Time associated with this file (high bytes). */
    uint32_t offset_lo;     /* Zero in the 1 frame/file case (low bytes) */
    uint32_t offset_hi;     /* Zero in the 1 frame/file case (high bytes) */
    uint32_t framesize_lo;  /* Number of bytes in frame (low bytes) */
    uint32_t framesize_hi;  /* Number of bytes in frame (high bytes) */
  } key_record_t;


  //// utility routines

  /*!
   * Extracts the low 32 bits of a 64 bit integer by masking.
   */
  uint32_t lobytes(const uint64_t& x) {
    uint32_t mask = 0xffffffff;
    return x & mask;
  }

  /*!
   * Extract the high 32 bits of a 64 bit integer by shifting.
   */
  uint32_t hibytes(const uint64_t& x) {
    return x >> 32;
  }

  /*!
   * Extract the low 32 bits of a 64 bit float as an integer.
   */
  uint32_t lobytes(const double& x) {
    union {
      uint64_t ival;
      double   dval;
    } u;
    u.dval = x;
    return lobytes(u.ival);
  }

  /*!
   * Extract the high 32 bits of a 64 bit float as an integer.
   */
  uint32_t hibytes(const double& x) {
    union {
      uint64_t ival;
      double   dval;
    } u;
    u.dval = x;
    return hibytes(u.ival);
  }

  /*!
   * The byte order associated with this machine.  We use
   * 1234 for little endian, 4321 for big endian, and
   * 3412 for the unlikely PDB endianism.
   */
  uint32_t machineEndianism() {
#if __BYTE_ORDER == __LITTLE_ENDIAN
    uint32_t byteorder = 1234;
#else
#if __BYTE_ORDER == __BIG_ENDIAN
    uint32_t byteorder = 4321;
#else
#ifdef PDB_ENDIAN
#if __BYTE_ORDER == __PDB_ENDIAN
    uint32_t byteorder = 3412;
#endif
#endif
#endif
#endif
    // If we get a compile error here, then __BYTE_ORDER
    // has an unexpected value.
    return byteorder;
  }

  uint64_t alignInteger( const uint64_t &x, unsigned border) {
    return x + (border - x%border)%border;
  }


  /*!
   * See RFC 1146 for Fletcher's Checksum (http://tools.ietf.org/html/rfc1146)
   */
  uint32_t fletcher( uint16_t *data, size_t len ) {
    uint32_t sum1 = 0xffff, sum2 = 0xffff;
 
    while (len) {
      unsigned tlen = len > 360 ? 360 : len;
      len -= tlen;
      do {
        sum1 += *data++;
        sum2 += sum1;
      } while (--tlen);
      sum1 = (sum1 & 0xffff) + (sum1 >> 16);
      sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    /* Second reduction step to reduce sums to 16 bits */
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    return sum2 << 16 | sum1;
  }

  /*!
   * Remove a file or directory.  For directories,
   * we recurse through subfiles and remove those
   * before attempting the ::rmdir();
   */
  void recursivelyRemove(std::string path) {
    struct stat statbuf;

    // -----------------------------------------------
    // Only try to unlink if the file exists
    // We recurse through directories and unlink
    // other files.
    // -----------------------------------------------

#ifdef DESRES_WIN32
    // Use ::stat instead of ::lstat on windows since there are no symlinks
    if (stat(path.c_str(),&statbuf) == 0) {
#else
    if (::lstat(path.c_str(),&statbuf) == 0) {
#endif
      if (!S_ISDIR(statbuf.st_mode)) {
        if (::unlink(path.c_str()) != 0) {
            throw std::runtime_error(strerror(errno));
        }
      } else {
        VMDDIR* directory = NULL;
        try {
          directory = vmd_opendir(path.c_str());
          if (directory) {
            // Remove subfiles
            char * entry;
            while( (entry=vmd_readdir(directory)) != NULL ) {
              // Don't unlink . or ..
              if (entry[0] == '.') {
                if (entry[1] == 0) continue;
                if (entry[1] == '.' && entry[2] == 0) continue;
              }
              recursivelyRemove(path + s_sep + entry);
            }
            vmd_closedir(directory);
            directory = NULL;

            // Remove the actual directory
            if (::rmdir(path.c_str()) != 0) {
              throw std::runtime_error(strerror(errno));
            }
          }
        } catch(...) {
          if (directory) vmd_closedir(directory);
          throw;
        }
      }
    }
  }


  ////////
  // CRC
  ////////
  
  typedef uint32_t crc;
  
  #define POLYNOMIAL 0x04C11DB7
  #define WIDTH  (8 * sizeof(crc))
  #define TOPBIT (1 << (WIDTH - 1))
  #define FINAL_XOR_VALUE 0xFFFFFFFF

  crc processByte( crc remainder, char msg ) {
        remainder ^= (msg << (WIDTH - 8));
        for (uint8_t bit = 8; bit > 0; --bit)
        {
            if (remainder & TOPBIT) {
                remainder = (remainder << 1) ^ POLYNOMIAL;
            } else {
                remainder = (remainder << 1);
            }
        }
        return remainder;
  }

  crc processBytes(const char *message, int nBytes) {
    crc  remainder = 0;	
    for (int byte = 0; byte < nBytes; ++byte) {
        remainder = processByte( remainder, message[byte] );
    }
    return remainder;
  } 

  int32_t cksum(const std::string &s) {
    ssize_t len = s.size();
    int32_t result = processBytes( s.c_str(), len );
  
    for ( ; len; len >>= 8) {
      result = processByte( result, len & 0xff );
    }
    return result ^ FINAL_XOR_VALUE;
  }

}


namespace {
  struct Blob {
    std::string type;
    uint64_t count;
    const void *data;
    bool byteswap;

    Blob() : count(0), data(0) {}
    Blob( const std::string &type_, uint64_t count_, const void *data_,
          uint32_t frame_endianism )
    : type(type_), count(count_), data(data_), byteswap(false) {
      uint32_t my_endianism = machineEndianism();
      if (frame_endianism != my_endianism) {
        if ( (frame_endianism==1234 && my_endianism==4321) || 
             (frame_endianism==4321 && my_endianism==1234) ) {
          byteswap=true;
        } else {
          throw std::runtime_error("Unable to handle frame endianness");
        }
      }
    }

    std::string str() const {
      if (type=="char" && count>0) {
        const char *s=(const char *)data;
        return std::string(s, s+count);
      }
      return "";
    }
    void get_float(float *buf) const {
      if (type=="float") {
        memcpy(buf, data, count*sizeof(float));
      } else if (type=="double") {
        const double *p = reinterpret_cast<const double *>(data);
        for (uint64_t i=0; i<count; i++) buf[i] = p[i];
      } else {
        memset(buf, 0, count*sizeof(float));
      }
      if (byteswap) swap4_unaligned(buf, count);
    }
    void get_double(double *buf) const {
      if (type=="double") {
        memcpy(buf, data, count*sizeof(double));
      } else if (type=="float") {
        const float *p = reinterpret_cast<const float *>(data);
        for (uint64_t i=0; i<count; i++) buf[i] = p[i];
      } else {
        memset(buf, 0, count*sizeof(double));
      }
      if (byteswap) swap8_unaligned(buf, count);
    }
    void get_int32(int32_t *buf) const {
      if (type=="int32_t") {
        memcpy(buf, data, count*sizeof(int32_t));
      } else {
        memset(buf, 0, count*sizeof(int32_t));
      }
      if (byteswap) swap4_unaligned(buf, count);
    }
    void get_uint32(uint32_t *buf) const {
      if (type=="uint32_t") {
        memcpy(buf, data, count*sizeof(uint32_t));
      } else {
        memset(buf, 0, count*sizeof(uint32_t));
      }
      if (byteswap) swap4_unaligned(buf, count);
    }
  };

  typedef std::map<std::string, Blob> BlobMap;
}

static const uint32_t magic_timekey = 0x4445534b;
static const uint32_t magic_frame   = 0x4445534d;
static const uint32_t s_version     = 0x00000100;
static const uint32_t s_irosetta    = 0x12345678;
static const float    s_frosetta    = 1234.5;
static const double   s_drosetta    = 1234.5e6;
static const uint32_t s_lrosetta_lo = 0x89abcdef;
static const uint32_t s_lrosetta_hi = 0x01234567;
static const uint32_t s_blocksize   = 4096;
static const uint32_t s_alignsize   = 8;

static uint64_t assemble64( uint32_t lo, uint32_t hi) {
  uint64_t hi64 = hi; 
  return (hi64 << 32) | lo; 
}

static double assembleDouble(uint32_t lo, uint32_t hi) {
  union {
    uint64_t ival;
    double   dval;
  } u;
  u.ival = assemble64(lo,hi);
  return u.dval;
}

static inline std::string addslash(const std::string& s){
    return (s.rbegin()[0] == '/') ? s : s + "/";
}

#define DD_RELPATH_MAXLEN (9) 
static std::string 
DDreldir(const std::string& fname, int ndir1, int ndir2){

    if( fname.find('/', 0) != std::string::npos ) {
      fprintf(stderr, "DDreldir: filename '%s' must not contain '/'\n",
          fname.c_str());
      return "";
    }

    uint32_t hash = cksum(fname);

    // uint32_t u1 = ndir1;
    // uint32_t u2 = ndir2;
    uint32_t d1, d2;
    char answer[DD_RELPATH_MAXLEN];
    if(ndir1 > 0){
	d1 = hash%ndir1;
	if(ndir2 > 0){
	    d2 = (hash/ndir1)%ndir2;
	    sprintf(answer, "%03x/%03x/", d1, d2);
	}else{
	    sprintf(answer, "%03x/", d1);
	}
    }else{
	sprintf(answer, "./");
    }
    return std::string(answer);
}

namespace {
  class DDException : public std::runtime_error{
  public:
      int eno;
      DDException(const std::string &text, int _eno=0) 
      : std::runtime_error(text + strerror(eno)), eno(_eno){}
  };
}

void DDmkdir(const std::string &dirpath, mode_t mode, int ndir1, int ndir2){
    std::string dpslash(addslash(dirpath));

    mode_t openmode = mode | 0300; // make sure we can write into the directory
    if( mkdir(dpslash.c_str(), openmode) < 0 )
	throw DDException("mkdir", errno);
	
    if( mkdir((dpslash + "not_hashed").c_str(), openmode) < 0 )
        throw DDException("mkdir not_hashed subdirectory", errno);

    FILE *fp = fopen((dpslash + "not_hashed/.ddparams").c_str(), "w");
    if(fp == NULL)
	throw DDException("fopen( .ddparams, \"w\" )", errno);
    if( fprintf(fp, "%d %d\n", ndir1, ndir2) < 0 ){
        fclose(fp);
	throw DDException("fprintf(.ddparams ...)", errno);
    }
    if( fclose(fp) )
	throw DDException("fclose(.ddparams)", errno);

    for(int i=0; i<ndir1; ++i){
	char sub[6];
	sprintf(sub, "%03x/", i);
	std::string dirsub = dpslash + sub;
        {
	    if( mkdir(dirsub.c_str(), openmode) < 0 )
		throw DDException("mkdir " + dirsub, errno);
	}
	for(int j=0; j<ndir2; ++j){
	    char subsub[6];
	    sprintf(subsub, "%03x", j);
	    std::string dirsubsub = dirsub + subsub;
	    if( mkdir(dirsubsub.c_str(), mode) < 0 ) // NOT openmode!
		throw DDException("mkdir " + dirsubsub, errno);
	}
        if( mode != openmode ){
            // change the mode back to what the user requested now
            // that we're done creating stuff...
            if( chmod(dirsub.c_str(), mode) < 0 )
                throw DDException("chmod " + dirsub, errno);
        }
    }
    if( mode != openmode ){
        // change the mode back to what the user requested now
        // that we're done creating stuff...
        if( chmod(dpslash.c_str(), mode) < 0 )
            throw DDException("chmod " + dpslash, errno);
        if( chmod((dpslash + "not_hashed").c_str(), mode) < 0 )
          throw DDException("chmod " + dpslash + "not_hashed", errno);
    }
}


static void 
DDgetparams(const std::string& dirpath, int *ndir1, int *ndir2) {
  // get ddparams, or assume (0,0) and let the frame file opens fail.
  *ndir1 = *ndir2 = 0;
  std::string dirslash(addslash(dirpath));
  // New convention - .ddparams is in not_hashed/.
  FILE *fp = fopen((dirslash + "not_hashed/.ddparams").c_str(), "r");
  // Allow the old convention of placing .ddparams in the top-level.
  if( fp == NULL && errno == ENOENT ) {
      fp = fopen((dirslash + ".ddparams").c_str(), "r");
  }
  if(fp != NULL) {
    if( fscanf(fp, "%d%d", ndir1, ndir2) != 2 ) 
      fprintf(stderr, "Failed to parse .ddparams; assuming flat structure\n");
    if( fclose(fp) ) {
      fprintf(stderr, "Warning: Failed to close .ddparams file: %s\n",
          strerror(errno));
    }
  }
}

static std::string framefile( const std::string &dtr,
                              size_t frameno, 
                              size_t frames_per_file,
                              int ndir1,
                              int ndir2) {
  unsigned frame_file = frameno / frames_per_file;
  std::ostringstream filename;
  filename << "frame" << std::setfill('0') << std::setw(9)
           << frame_file;
  std::string fname = filename.str();

  std::string fullpath(dtr);
  fullpath += "/";
  fullpath += DDreldir(fname, ndir1, ndir2);
  fullpath += fname;
  return fullpath;
}

static BlobMap read_frame( const void *mapping ) {

    const char *base = reinterpret_cast<const char *>(mapping);
    const header_t *header = reinterpret_cast<const header_t*>(base);
    uint32_t size_header_block = ntohl(header->size_header_block);
    uint32_t frames_endianism = ntohl(header->endianism);
    uint32_t frames_nlabels = ntohl(header->nlabels);
    uint32_t size_meta_block = ntohl(header->size_meta_block);
    uint32_t size_typename_block = ntohl(header->size_typename_block);
    uint32_t size_label_block = ntohl(header->size_label_block);
    uint32_t size_scalar_block = ntohl(header->size_scalar_block);
    //uint32_t size_field_block_lo = ntohl(header->size_field_block_lo);
    //uint32_t size_field_block_hi = ntohl(header->size_field_block_hi);
    //uint64_t size_field_block = assemble64(size_field_block_lo,
                                           //size_field_block_hi);

    uint64_t offset_header_block = 0;
    uint64_t offset_meta_block = offset_header_block + size_header_block;
    uint64_t offset_typename_block = offset_meta_block + size_meta_block;
    uint64_t offset_label_block = offset_typename_block + size_typename_block;
    uint64_t offset_scalar_block = offset_label_block + size_label_block;
    uint64_t offset_field_block = offset_scalar_block + size_scalar_block;
    //uint64_t offset_crc_block = offset_field_block + size_field_block;

    const metadisk_t* diskmeta  = reinterpret_cast<const metadisk_t*>(base+offset_meta_block);
    const char* typenames = reinterpret_cast<const char*>(base+offset_typename_block);
    const char* labels    = reinterpret_cast<const char*>(base+offset_label_block); 
    const char* scalars   = reinterpret_cast<const char*>(base+offset_scalar_block);
    const char* fields    = reinterpret_cast<const char*>(base+offset_field_block);

    std::vector<std::string> types;
    while(*typenames) {
      if (typenames >= labels) {
        fprintf(stderr, "More typenames than labels!\n");
        break;
      }
      std::string type(typenames);
      types.push_back(type);
      typenames += type.size()+1;
    }

    BlobMap blobs;

    for (size_t ii=0; ii<frames_nlabels; ++ii) {
      std::string label(labels);
      labels += label.size()+1;
      // Pull out the typecode, elementsize, and count
      uint32_t code = ntohl(diskmeta[ii].type);
      uint32_t elementsize = ntohl(diskmeta[ii].elementsize);
      uint32_t count_lo = ntohl(diskmeta[ii].count_lo);
      uint32_t count_hi = ntohl(diskmeta[ii].count_hi);
      uint64_t count = assemble64(count_lo,count_hi);
      uint64_t nbytes = elementsize*count;

      const char *addr=0;
      if (count <= 1) {
        addr=scalars;
        scalars += alignInteger(nbytes, s_alignsize);
      } else {
        addr=fields;
        fields += alignInteger(nbytes, s_alignsize);
      }
      try {
        blobs[label] = Blob( types.at(code), count, addr, frames_endianism );
      }
      catch (std::exception &e) {
        fprintf(stderr, "Failed fetching '%s' data from frame\n", 
            label.c_str());
      }
    }
    return blobs;
}

static void *map_file( int fd, off_t offset, size_t &framesize ) {
  if (fd<=0) {
    fprintf(stderr, "map_file: bad file descriptor\n");
    return MAP_FAILED;
  }
  if (framesize==0) {
    struct stat statbuf;
    if (fstat(fd,&statbuf)!=0) {
      fprintf(stderr, "Could not stat file: %s\n", strerror(errno));
      return MAP_FAILED;
    }
    framesize=statbuf.st_size-offset;
  }

  void *mapping = ::mmap(0, framesize, PROT_READ|PROT_WRITE,MAP_PRIVATE,
                         fd, offset);
  if (mapping==MAP_FAILED) {
    fprintf(stderr, "Failed to map file: %s\n", strerror(errno));
  }
  return mapping;
}

namespace {

}

static void get_rmass(const std::string &metafile, std::vector<float> &rmass) {

  rmass.resize(0);
  int meta_fd = open(metafile.c_str(), O_RDONLY|O_BINARY);
  size_t framesize=0;
  void *meta_mapping = map_file( meta_fd, 0, framesize );
  if (meta_mapping==MAP_FAILED) {
    close(meta_fd);
    return;
  }
  BlobMap meta_blobs = read_frame( meta_mapping );
  if (meta_blobs.find("INVMASS")!=meta_blobs.end()) {

    Blob blob=meta_blobs["INVMASS"];
    rmass.resize(blob.count);
    blob.get_float(&rmass[0]);

    // have to permute if GID is present
    if (meta_blobs.find("GID")!=meta_blobs.end()) {
      blob=meta_blobs["GID"];
      std::vector<uint32_t> gids;
      gids.resize(blob.count);
      blob.get_uint32(&gids[0]);
    
      std::vector<float> perm_rmass(rmass.size());
      // assuming gids are contiguous, 0-based.
      for (unsigned i=0; i<gids.size(); i++) {
        perm_rmass[gids[i]] = rmass[i];
      }
      rmass = perm_rmass;
    }
  }
  ::munmap(MMAPTYPE(meta_mapping), framesize);
  close(meta_fd);
}


///// Plugin interface
#include <molfile_plugin.h>

namespace {
  struct Handle {
    std::string dtr;
    std::vector<key_record_t> keys;  // one per frame
    std::vector<float> rmass;   // could be empty vector!
    const uint32_t natoms;
    const uint32_t frames_per_file;
    const int ndir1;
    const int ndir2;
    uint32_t curframe;
    bool with_velocity;

    int frame_fd; // for writing

    Handle( const std::string &path, uint32_t natoms_,
            uint32_t frames_per_file_, int n1, int n2 )
    : dtr(path), natoms(natoms_), frames_per_file(frames_per_file_),
      ndir1(n1), ndir2(n2),
      curframe(0), with_velocity(false)
    {}

    ~Handle() { }
  };
}

static void *
open_file_read( const char *filename, const char *filetype, int *natoms ) {

  // check for "clickme.dtr"
  std::string fname(filename);
  std::string::size_type pos = fname.rfind( "clickme.dtr" );
  if (pos != std::string::npos) {
    fname.resize(pos);
    filename = fname.c_str();
  }

  /* Read the timekeys file */
  const std::string dtr = filename;
  FILE *fd=0;
  std::string timekeys_path = dtr;
  timekeys_path += s_sep;
  timekeys_path += "timekeys";
  fd = fopen( timekeys_path.c_str(), "rb" );
  if (!fd) {
    fprintf(stderr, "Could not find timekeys file at %s\n", 
        timekeys_path.c_str());
    return NULL;
  }


  /* check the magic number */
  key_prologue_t prologue[1];
  if (fread( prologue, sizeof(key_prologue_t), 1, fd )!=1) {
    fprintf(stderr, "Failed to read key prologue from %s\n",
        timekeys_path.c_str());
    fclose(fd);
    return NULL;
  }
  prologue->magic = htonl(prologue->magic);
  if (prologue->magic != magic_timekey) {
    fprintf(stderr, "timekeys magic number %x doesn't match %x\n",
        prologue->magic, magic_timekey);
    fclose(fd);
    return NULL;
  }

  /* get frames per file and key record size */
  prologue->frames_per_file = ntohl( prologue->frames_per_file );
  prologue->key_record_size = ntohl( prologue->key_record_size );

  /* read all key records */
  fseek(fd, 0, SEEK_END);
  off_t keyfile_size = ftello(fd);
  size_t nframes = (keyfile_size-sizeof(key_prologue_t))/sizeof(key_record_t);

  if (!nframes) {
    fprintf(stderr, "Error, empty trajectory\n");
    fclose(fd);
    return NULL;
  }

  std::vector<key_record_t> keys(nframes);
  fseek(fd, sizeof(key_prologue_t), SEEK_SET);
  if (fread(&keys[0], sizeof(key_record_t), nframes, fd)!=nframes) {
    fprintf(stderr, "Failed to read all timekeys records: %s\n",
        strerror(errno));
    fclose(fd);
    return NULL;
  }
  fclose(fd);

  int n1, n2;
  DDgetparams(dtr, &n1, &n2);

  // read the first frame to see how many atoms there are, and whether 
  // there are any velocities.
  *natoms = -1;
  bool has_velocities=false;
  {
    std::string fname=framefile(dtr, 0, prologue->frames_per_file, n1, n2);
    int fd = open(fname.c_str(), O_RDONLY|O_BINARY);
    size_t framesize=0;
    unsigned i;
    void *mapping = map_file( fd, 0, framesize );
    if (mapping==MAP_FAILED)  {
      fprintf(stderr, "Failed to find frame at %s\n", fname.c_str());
      close(fd);
      return NULL;
    }
    BlobMap blobs = read_frame(mapping);

    // I'm aware of three possible sources of positions: 
    //  "POSN" (the original frameset format)
    //  "POSITION" (the wrapped frameset formats)
    //  "POS" (anton trajectories)
    const char *posnames[] = { "POSN", "POSITION", "POS" };
    for (i=0; i<3; i++) {
      if (blobs.find(posnames[i])!=blobs.end()) {
        *natoms = blobs[posnames[i]].count / 3;
        break;
      }
    }
    // similar for velocities
    const char *velnames[] = { "MOMENTUM", "VELOCITY" };
    for (i=0; i<2; i++) {
      if (blobs.find(velnames[i])!=blobs.end()) {
        has_velocities=true;
        break;
      }
    }
    ::munmap(MMAPTYPE(mapping), framesize);
    close(fd);
  }


  Handle *handle = new Handle(dtr, *natoms, prologue->frames_per_file, n1, n2);
  handle->keys = keys;

  // get the reciprocal mass, if present
  get_rmass( dtr + s_sep + "metadata", handle->rmass );
  handle->with_velocity = has_velocities && !strcmp(filetype, "dtrv");
  return handle;
}

static int 
read_timestep_metadata(void *v, molfile_timestep_metadata *m) {
  Handle* handle = reinterpret_cast<Handle*>(v);
  m->has_velocities = handle->with_velocity;
  m->count = handle->keys.size();
  return MOLFILE_SUCCESS;
}

static double dotprod(const double *x, const double *y) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

static void read_homebox( const double *box,
                          molfile_timestep_t *ts ) {

  ts->A = ts->B = ts->C = 1;
  ts->alpha = ts->beta = ts->gamma = 90;

  double A[3] = { box[0], box[3], box[6] };
  double B[3] = { box[1], box[4], box[7] };
  double C[3] = { box[2], box[5], box[8] };

  // store lengths
  ts->A = sqrt(dotprod(A,A));
  ts->B = sqrt(dotprod(B,B));
  ts->C = sqrt(dotprod(C,C));

  // compute angles
  double cosAB = dotprod(A,B)/(ts->A * ts->B);
  double cosAC = dotprod(A,C)/(ts->A * ts->C);
  double cosBC = dotprod(B,C)/(ts->B * ts->C);

  // clamp
  if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
  if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
  if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;

  // convert to angles using asin to avoid nasty rounding when we are
  // close to 90 degree angles.
  ts->alpha = 90.0 - asin(cosBC) * 90.0 / M_PI_2; /* cosBC */
  ts->beta  = 90.0 - asin(cosAC) * 90.0 / M_PI_2; /* cosAC */
  ts->gamma = 90.0 - asin(cosAB) * 90.0 / M_PI_2; /* cosAB */
}

void write_homebox( const molfile_timestep_t * ts,
                    float * box ) {

  double A[3], B[3], C[3];

  // Convert VMD's unit cell information
  double cosBC = sin( ((90 - ts->alpha ) / 180) * M_PI );
  double cosAC = sin( ((90 - ts->beta  ) / 180) * M_PI );
  double cosAB = sin( ((90 - ts->gamma ) / 180) * M_PI );
  double sinAB = cos( ((90 - ts->gamma ) / 180) * M_PI );

  double Ax = ts->A;
  double Ay = 0;
  double Az = 0;
  double Bx = ts->B * cosAB;
  double By = ts->B * sinAB;
  double Bz = 0;
  double Cx,Cy,Cz;
  if (sinAB != 0) {
    Cx = cosAC;
    Cy = (cosBC - cosAC*cosAB) / sinAB;
    Cz = sqrt(1-Cx*Cx-Cy*Cy);
    Cx *= ts->C;
    Cy *= ts->C;
    Cz *= ts->C;
  } else {
    Cx=Cy=Cz=0;
  }
  A[0] = Ax; A[1] = Ay; A[2] = Az;
  B[0] = Bx; B[1] = By; B[2] = Bz;
  C[0] = Cx; C[1] = Cy; C[2] = Cz;

  // put vectors in column of homebox
  box[0] = A[0]; box[3] = A[1]; box[6] = A[2];
  box[1] = B[0]; box[4] = B[1]; box[7] = B[2];
  box[2] = C[0]; box[5] = C[1]; box[8] = C[2];
}

static int handle_wrapped_v2(const Handle *h,
                             BlobMap &blobs,
                             molfile_timestep_t *ts) {

  // just read POSITION in either single or double precision
  if (blobs.find("POSITION")==blobs.end()) {
    fprintf(stderr, "ERROR, Missing POSITION field in frame %d\n", h->curframe);
    return MOLFILE_ERROR;
  }
  Blob pos=blobs["POSITION"];
  if (pos.count != 3*h->natoms) {
    fprintf(stderr, "ERROR, Expected %d elements in POSITION; got %ld\n",
        3*h->natoms, pos.count);
    return MOLFILE_ERROR;
  }
  pos.get_float(ts->coords);

  if (h->with_velocity && ts->velocities && blobs.find("VELOCITY")!=blobs.end()) {
    Blob vel=blobs["VELOCITY"];
    if (vel.count != 3*h->natoms) {
      fprintf(stderr, "ERROR, Expected %d elements in VELOCITY; got %ld\n",
          3*h->natoms, vel.count);
      return MOLFILE_ERROR;
    }
    vel.get_float(ts->velocities);
  }

  if (blobs.find("UNITCELL")!=blobs.end()) {
    double box[9];
    blobs["UNITCELL"].get_double(box);
    read_homebox( box, ts );
  }
  return MOLFILE_SUCCESS;
}

namespace {

  inline void
  compute_center(int partition,
                 int nx, int ny, int nz,
                 float b0, float b1, float b2,
                 float b3, float b4, float b5,
                 float b6, float b7, float b8,
                 float* cx, float* cy, float* cz) {
    double nu, nv, nw, mu, mv, mw;
    double xc, yc, zc;

    // -----------------------------------------------
    // Map the partition number to its "mesh" position
    // (see define_mesh_collective in topology.c)
    // -----------------------------------------------
    int hmx = partition;
    int hmy  = hmx / nx;     /* y = y + ny*( z + nz*r ) */
    int hmz  = hmy / ny;     /* z = z + nz*r */
    hmx -= hmy * nx;         /* x = x */
    hmy -= hmz * ny;         /* y = y */

    nu = (double)nx;
    nv = (double)ny;
    nw = (double)nz;

    // -----------------------------------------------
    // Code adapted from configure_global_cell in
    // topology.c
    // -----------------------------------------------
    mu = -0.5*(nu-1) + (double)hmx;
    mv = -0.5*(nv-1) + (double)hmy;
    mw = -0.5*(nw-1) + (double)hmz;

    // We used to do FORCE_PRECISION(xc,float) here, but that
    // seems unnecessary in the context of trajectory writing.
    xc = b0*mu + b1*mv + b2*mw; 
    yc = b3*mu + b4*mv + b5*mw; 
    zc = b6*mu + b7*mv + b8*mw; 

    *cx = xc;
    *cy = yc;
    *cz = zc;
  }

  inline int 
  posn_momentum_v_1(int32_t nx, int32_t ny, int32_t nz,
                    uint64_t nparticles,
                    const double  * home_box,
                    const uint32_t* gid,
                    const uint32_t* npp,
                    const float   * rm, // reciprocal mass
                    const float* posn, const float* momentum,
                    /* returns */
                    float *position, float *velocity, double *box) {

    // bounding box is a straight multiple of the home box
    if (box) {
      box[0] = home_box[0]*nx;
      box[1] = home_box[1]*ny;
      box[2] = home_box[2]*nz;
        
      box[3] = home_box[3]*nx;
      box[4] = home_box[4]*ny;
      box[5] = home_box[5]*nz;

      box[6] = home_box[6]*nx;
      box[7] = home_box[7]*ny;
      box[8] = home_box[8]*nz;
    }


    int partition = 0;
    int remaining = 0;
    float cx = 0;
    float cy = 0;
    float cz = 0;
    float ux = home_box[0];
    float vx = home_box[1];
    float wx = home_box[2];
    float uy = home_box[3];
    float vy = home_box[4];
    float wy = home_box[5];
    float uz = home_box[6];
    float vz = home_box[7];
    float wz = home_box[8];

    for(uint64_t i=0; i<nparticles; ++i) {
      if (remaining == 0) {
        do {
          remaining = npp[partition];
          ++partition;
        } while (!remaining); // skip empty partitions
          compute_center(partition-1, nx,ny,nz, ux,vx,wx, uy,vy,wy, uz,vz,wz,
                        &cx,&cy,&cz);
      }
      uint32_t id = gid[i];
      if (id >= nparticles) {
        fprintf(stderr, "non-contiguous particles\n");
        return MOLFILE_ERROR;
      }

      if (posn) {
        float x = posn[3*i+0];
        float y = posn[3*i+1];
        float z = posn[3*i+2];

        position[3*id+0] = ux*x + vx*y + wx*z + cx;
        position[3*id+1] = uy*x + vy*y + wy*z + cy;
        position[3*id+2] = uz*x + vz*y + wz*z + cz;
      }

      if (velocity && momentum && rm) {
        velocity[3*id+0] = momentum[3*i+0]*rm[id];
        velocity[3*id+1] = momentum[3*i+1]*rm[id];
        velocity[3*id+2] = momentum[3*i+2]*rm[id];
      } else if (velocity) {
        velocity[3*id+0] = 0.0;
        velocity[3*id+1] = 0.0;
        velocity[3*id+2] = 0.0;
      }
      --remaining;
    }
    return MOLFILE_SUCCESS;
  }
}

static int handle_posn_momentum_v1(const Handle *h,
                                   BlobMap &blobs,
                                   molfile_timestep_t *ts) {

  int32_t nx, ny, nz;
  double home_box[9], box[9];
  blobs["HOME_BOX"].get_double(home_box);
  blobs["NX"].get_int32(&nx);
  blobs["NY"].get_int32(&ny);
  blobs["NZ"].get_int32(&nz);
  
  std::vector<uint32_t> gid, npp;
  std::vector<float> pos, mtm;
  Blob gidblob=blobs["GID"];
  Blob nppblob=blobs["NPP"];
  Blob posblob=blobs["POSN"];
  Blob mtmblob=blobs["MOMENTUM"];

  if (gidblob.count != h->natoms) {
    fprintf(stderr, "Missing GID field\n");
    return MOLFILE_ERROR;
  }
  if (posblob.count != 3*h->natoms) {
    fprintf(stderr, "Missing POSN field\n");
    return MOLFILE_ERROR;
  }
  gid.resize(gidblob.count);
  npp.resize(nppblob.count);
  pos.resize(posblob.count);
  mtm.resize(mtmblob.count);

  gidblob.get_uint32(&gid[0]);
  nppblob.get_uint32(&npp[0]);
  posblob.get_float(&pos[0]);

  if (h->rmass.size() && h->with_velocity) mtmblob.get_float(&mtm[0]);

  posn_momentum_v_1( nx, ny, nz, h->natoms, home_box, 
                     &gid[0], &npp[0], 
                     h->rmass.size() ? &h->rmass[0] : NULL,
                     &pos[0],
                     &mtm[0],
                     ts->coords,
                     ts->velocities,
                     box );

  read_homebox( box, ts );
  return MOLFILE_SUCCESS;
}

static int handle_wrapped_v1(const Handle *h,
                             BlobMap &blobs,
                             molfile_timestep_t *ts) {
  {
    // homebox
    double home_box[9], box[9];
    int32_t nx, ny, nz;
    blobs["HOME_BOX"].get_double(home_box);
    blobs["NX"].get_int32(&nx);
    blobs["NY"].get_int32(&ny);
    blobs["NZ"].get_int32(&nz);
    box[0] = home_box[0]*nx;
    box[1] = home_box[1]*ny;
    box[2] = home_box[2]*nz;
      
    box[3] = home_box[3]*nx;
    box[4] = home_box[4]*ny;
    box[5] = home_box[5]*nz;
 
    box[6] = home_box[6]*nx;
    box[7] = home_box[7]*ny;
    box[8] = home_box[8]*nz;
    read_homebox( box, ts );
  }

  Blob posblob=blobs["POSN"];
  Blob velblob=blobs["VELOCITY"];

  // get positions
  if (posblob.count != 3*h->natoms) {
    fprintf(stderr, "Missing POSN field\n");
    return MOLFILE_ERROR;
  }
  posblob.get_float(ts->coords);
  
  // if required, get velocities
  if (ts->velocities && velblob.count > 0) {
    if (velblob.count != 3*h->natoms) {
      fprintf(stderr, "VELOCITY field has %d values; expected %d\n",
          velblob.count, 3*h->natoms);
      return MOLFILE_ERROR;
    }
    velblob.get_float(ts->velocities);
  }
  return MOLFILE_SUCCESS;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  Handle *h = reinterpret_cast<Handle *>(v);
  if (h->curframe >= h->keys.size()) return MOLFILE_EOF;
  if (!ts) {
    ++h->curframe;
    return MOLFILE_SUCCESS;
  }
  uint32_t iframe = h->curframe;
  const std::vector<key_record_t> &keys = h->keys;
  int rc = MOLFILE_SUCCESS;
  {
    off_t offset=0;
    size_t framesize=0;
    if (h->frames_per_file != 1) {
      offset = assemble64( ntohl(keys[iframe].offset_lo), 
                           ntohl(keys[iframe].offset_hi) );
      framesize = assemble64( ntohl(keys[iframe].framesize_lo), 
                              ntohl(keys[iframe].framesize_hi) );

    }
    ts->physical_time = assembleDouble(ntohl(keys[iframe].time_lo),
                                       ntohl(keys[iframe].time_hi));

    std::string fname = framefile( h->dtr, iframe, h->frames_per_file, 
        h->ndir1, h->ndir2 );
    int fd = open(fname.c_str(), O_RDONLY|O_BINARY);
    void *mapping = map_file( fd, offset, framesize );
    BlobMap blobs = read_frame(mapping);

    // Now, dispatch to routines based on format
    std::string format = blobs["FORMAT"].str();
    if (format=="WRAPPED_V_2" || format == "DBL_WRAPPED_V_2") {
      rc = handle_wrapped_v2(h, blobs, ts);

    } else if (format=="POSN_MOMENTUM_V_1" || format=="DBL_POSN_MOMENTUM_V_1") {
      rc = handle_posn_momentum_v1(h, blobs, ts);

    } else if (format=="WRAPPED_V_1" || format == "DBL_WRAPPED_V_1") {
      rc = handle_wrapped_v1(h, blobs, ts);

    } else {
      fprintf(stderr, "ERROR, can't handle format %s\n", format.c_str());
      rc = MOLFILE_ERROR;
    }
    munmap(MMAPTYPE(mapping), framesize);
    close(fd);
  }
  ++h->curframe;
  return rc;
}


static void close_file_read( void *v ) {
  Handle *h = reinterpret_cast<Handle *>(v);
  delete h;
}

static void *open_file_write(const char *path, const char *type, int natoms) {

  Handle *h = NULL;
  try {
    std::string m_directory(path);
    char cwd[4096];

    while(m_directory.size() > 0 && m_directory[m_directory.size()-1] == s_sep) {
      m_directory.erase(m_directory.size()-1);
    }

    if ( m_directory[0] != s_sep) {
      if (! ::getcwd(cwd,sizeof(cwd))) {
        throw std::runtime_error(strerror(errno));
      }
      m_directory = std::string(cwd) + s_sep + m_directory;
    }

    recursivelyRemove(m_directory);
    ::DDmkdir(m_directory,0777,0, 0);

    // craft an empty metadata frame
    std::string metadata_file = m_directory + s_sep + "metadata";
    FILE *fd = fopen(metadata_file.c_str(), "wb");

    const uint64_t framesize = 4096;
    const uint64_t size_header_block = 
      alignInteger( sizeof(header_t), s_alignsize );

    std::vector<char> bytes(framesize);
    char *base = &bytes[0];
    memset( base, 0, framesize );

    header_t *header = reinterpret_cast<header_t*>(base);
    memset(header,0,sizeof(header_t));
    header->required.magic = htonl(magic_frame);
    header->required.version = htonl(s_version);

    header->required.framesize_lo = htonl(lobytes(framesize));
    header->required.framesize_hi = htonl(hibytes(framesize));

    header->size_header_block = htonl(size_header_block);
    header->unused0 = 0;
    uint64_t lrosetta = assemble64(s_lrosetta_lo,s_lrosetta_hi);
    header->irosetta = s_irosetta;
    header->frosetta = s_frosetta;

    header->drosetta_lo = lobytes(s_drosetta);
    header->drosetta_hi = hibytes(s_drosetta);
    header->lrosetta_lo = lobytes(lrosetta);
    header->lrosetta_hi = hibytes(lrosetta);

    header->endianism = htonl(machineEndianism());

    fwrite( base, framesize, 1, fd );
    fclose(fd);

    h = new Handle( m_directory, natoms, 1, 0, 0 );
    h->with_velocity = !strcmp(type, "dtrv");

    // we're going to write all timesteps to a single frame
    std::string filepath = framefile( m_directory, 0, 1, 0, 0 );
    h->frame_fd = open(filepath.c_str(),O_WRONLY|O_CREAT|O_BINARY,0666);
    if (h->frame_fd<0) throw std::runtime_error(strerror(errno));
  }
  catch (std::exception &e) {
    delete h; h=NULL;
    fprintf(stderr, "%s\n", e.what());
    return NULL;
  }
  return h;
}

namespace {
  struct meta_t {
    std::string label;
    std::string typecode;
    uint32_t elementsize;
    uint64_t count;
    const char *bytes;
    meta_t() {}
    meta_t(const std::string &l, const std::string &t, uint32_t e, uint32_t c,
           const void *b)
    : label(l), typecode(t), elementsize(e), count(c), 
    bytes(reinterpret_cast<const char *>(b)) {}
  };
  typedef std::vector<meta_t> MetaList;

  uint64_t typename_size(const MetaList &meta) {
    // just the set of distinct types
    uint64_t sz=0;
    typedef std::set<std::string> Typemap;
    Typemap types;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      types.insert(m->typecode);
    for (Typemap::const_iterator s=types.begin(); s!=types.end();++s)
      sz += s->size() + 1;
    sz += 1;
    return alignInteger(sz, s_alignsize);
  }

  uint64_t label_size(const MetaList &meta) {
    uint64_t sz=0;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      sz += m->label.size() + 1;
    sz += 1;
    return alignInteger(sz, s_alignsize);
  }

  uint64_t scalar_size(const MetaList &meta) {
    uint64_t sz=0;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      if (m->count <= 1) 
        sz += alignInteger( m->elementsize * m->count, s_alignsize );
    return sz;
  }
  uint64_t field_size(const MetaList &meta) {
    uint64_t sz=0;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      if (m->count > 1) 
        sz += alignInteger( m->elementsize * m->count, s_alignsize );
    return sz;
  }
}

static int write_timestep(void *v, const molfile_timestep_t *ts) {
  Handle *h = reinterpret_cast<Handle *>(v);
  try {

    int natoms = h->natoms;
    bool with_velocity = h->with_velocity;

    static const char *format = "WRAPPED_V_2";
    static const char *title = "written by VMD";
    float box[9];
    write_homebox( ts, box );

    double time = ts->physical_time;

    std::vector<meta_t> meta;
    meta.push_back( 
        meta_t( "FORMAT", "char", sizeof(char), strlen(format), format ));
    meta.push_back( 
        meta_t( "TITLE", "char", sizeof(char), strlen(title), title ));
    meta.push_back( 
        meta_t( "CHEMICAL_TIME", "double", sizeof(double), 1, &time));
    meta.push_back( 
        meta_t( "UNITCELL", "float", sizeof(float), 9, box ));
    meta.push_back( 
        meta_t( "POSITION", "float", sizeof(float), 3*natoms, ts->coords ));
    if (with_velocity) meta.push_back( 
        meta_t( "VELOCITY", "float", sizeof(float), 3*natoms, ts->velocities ));

    uint64_t offset_header_block = 0;
    uint64_t size_header_block =
      alignInteger( sizeof(header_t), s_alignsize );

    uint64_t offset_meta_block = offset_header_block + size_header_block;
    uint64_t size_meta_block = 
      alignInteger( meta.size()*sizeof(metadisk_t), s_alignsize );

    uint64_t offset_typename_block = offset_meta_block + size_meta_block;
    uint64_t size_typename_block = typename_size(meta);

    uint64_t offset_label_block = offset_typename_block + size_typename_block;
    uint64_t size_label_block = label_size(meta);

    uint64_t offset_scalar_block = offset_label_block + size_label_block;
    uint64_t size_scalar_block = scalar_size(meta);

    uint64_t offset_field_block = offset_scalar_block + size_scalar_block;
    uint64_t size_field_block = field_size(meta);

    uint64_t offset_crc_block = offset_field_block + size_field_block;
    uint64_t size_crc_block = sizeof(uint32_t);

    uint64_t offset_padding_block = offset_crc_block + size_crc_block;
    uint64_t size_padding_block = 
      alignInteger(offset_padding_block,s_blocksize) - offset_padding_block;

    uint64_t framesize = offset_padding_block + size_padding_block;

    // construct the frame
    char *base = new char[framesize];
    memset( base, 0, framesize );

    header_t *header = reinterpret_cast<header_t*>(base+offset_header_block);
    metadisk_t* diskmeta  = reinterpret_cast<metadisk_t*>(base+offset_meta_block);
    char*       typenames = reinterpret_cast<char*>(base+offset_typename_block);
    char*       labels    = reinterpret_cast<char*>(base+offset_label_block);
    char*       scalars   = reinterpret_cast<char*>(base+offset_scalar_block);
    char*       fields    = reinterpret_cast<char*>(base+offset_field_block);
    uint32_t*   crc       = reinterpret_cast<uint32_t*>(base+offset_crc_block);
    //char*       padding   = reinterpret_cast<char*>(base+offset_padding_block);

    /*** header ***/
    memset(header,0,sizeof(header_t));
    header->required.magic = htonl(magic_frame);
    header->required.version = htonl(s_version);

    header->required.framesize_lo = htonl(lobytes(framesize));
    header->required.framesize_hi = htonl(hibytes(framesize));

    header->size_header_block = htonl(size_header_block);
    header->unused0 = 0;
    uint64_t lrosetta = assemble64(s_lrosetta_lo,s_lrosetta_hi);
    header->irosetta = s_irosetta;
    header->frosetta = s_frosetta;

    header->drosetta_lo = lobytes(s_drosetta);
    header->drosetta_hi = hibytes(s_drosetta);
    header->lrosetta_lo = lobytes(lrosetta);
    header->lrosetta_hi = hibytes(lrosetta);

    header->endianism = htonl(machineEndianism());
    header->nlabels = htonl(meta.size());
    header->size_meta_block = htonl(size_meta_block);
    header->size_typename_block = htonl(size_typename_block);

    header->size_label_block = htonl(size_label_block);
    header->size_scalar_block = htonl(size_scalar_block);
    header->size_field_block_lo = htonl(lobytes(size_field_block));
    header->size_field_block_hi = htonl(hibytes(size_field_block));

    header->size_crc_block = htonl(size_crc_block);
    header->size_padding_block = htonl(size_padding_block);
    header->unused1 = 0;
    header->unused2 = 0;

    std::map<std::string,unsigned> typemap;

    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m) {

      if (typemap.find(m->typecode)==typemap.end()) {
        unsigned code=typemap.size();
        typemap[m->typecode]=code;
        typenames=std::copy(m->typecode.begin(), m->typecode.end(), typenames);
        *typenames++ = 0;
      }

      diskmeta->type = htonl( typemap[m->typecode] );
      diskmeta->elementsize = htonl( m->elementsize );
      diskmeta->count_lo = htonl( lobytes( m->count ));
      diskmeta->count_hi = htonl( hibytes( m->count ));
      diskmeta++;

      labels=std::copy(m->label.begin(), m->label.end(), labels);
      *labels++ = 0;

      uint64_t nbytes = m->count*m->elementsize;
      if (m->count <= 1) {
        memcpy( scalars, m->bytes, nbytes );
        scalars += alignInteger( nbytes, s_alignsize );
      } else {
        memcpy( fields, m->bytes, nbytes );
        fields += alignInteger( nbytes, s_alignsize );
      }
    }
    *crc = fletcher(reinterpret_cast<uint16_t*>(base),offset_crc_block/2);
    

    // write the data to disk
    ssize_t n;
    do {
      n = ::write(h->frame_fd, base, framesize );
    } while (n<0 && errno == EINTR);
    delete [] base;
    if (n<0)
      throw std::runtime_error(std::string("writing frame: ")+strerror(errno));

    uint64_t framefile_offset = h->keys.size() * framesize;

    // add an entry to the keyfile list
    key_record_t timekey;
    timekey.time_lo = htonl(lobytes(time));
    timekey.time_hi = htonl(hibytes(time));
    timekey.offset_lo = htonl(lobytes(framefile_offset));
    timekey.offset_hi = htonl(hibytes(framefile_offset));
    timekey.framesize_lo = htonl(lobytes(framesize));
    timekey.framesize_hi = htonl(hibytes(framesize));
    h->keys.push_back(timekey);
  } 
  catch (std::exception &e) {
    fprintf(stderr, "Write failed: %s\n",e.what());
    return MOLFILE_ERROR;
  }
  return MOLFILE_SUCCESS;
}

static void close_file_write( void * v ) {
  Handle *h = reinterpret_cast<Handle *>(v);

  std::string timekeys_path = h->dtr + s_sep + "timekeys";
  FILE *fd = fopen( timekeys_path.c_str(), "wb" );
  if (!fd) {
    fprintf(stderr, "Opening timekeys failed: %s\n", strerror(errno));
  } else {
    key_prologue_t prologue[1];
    prologue->magic = htonl(magic_timekey);
    prologue->frames_per_file = htonl(h->keys.size());
    prologue->key_record_size = htonl(sizeof(key_record_t));
    fwrite( prologue, sizeof(key_prologue_t), 1, fd );
    fwrite( &h->keys[0], sizeof(key_record_t), h->keys.size(), fd );
    fclose(fd);
  }
  ::close(h->frame_fd);
  delete h;
}

#if defined(TEST_DTRPLUGIN)

int main(int argc, char *argv[]) {

  /* check input arguments */
  if (argc != 3) {
    fprintf(stderr, "Usage: %s input.dtr output.dtr\n", argv[0]);
    return 1;
  }
  int natoms;
  
  /* read in all frames */
  void *handle = open_file_read( argv[1], "dtr", &natoms);
  printf("got %d atoms\n", natoms);
  molfile_timestep_t ts[1];
  std::vector<float> timesteps;
  do {
    timesteps.resize( timesteps.size()+3*natoms );
    ts->coords = &timesteps[ timesteps.size() - 3*natoms ];
    ts->velocities = NULL;
  } while (read_next_timestep(handle, natoms, ts)==MOLFILE_SUCCESS);
  int nframes = timesteps.size()/(3*natoms) - 1;

  printf("read %d steps\n", nframes );
  close_file_read(handle);

  /* write out all frames */
  handle = open_file_write( argv[2], "dtr", natoms);
  if (!handle) return 1;
  ts->coords = &timesteps[0];
  for (int i=0; i<nframes; i++) {
    if (write_timestep( handle, ts )!=MOLFILE_SUCCESS) {
      fprintf(stderr, "failed to write timestep %d/%d\n", i+1, nframes);
      return 1;
    }
    ts->coords += 3*natoms;
  }
  printf("wrote %d steps\n", nframes);
  close_file_write(handle);

  /* now try to read it back in */
  int new_natoms;
  handle = open_file_read( argv[2], "dtr", &new_natoms);
  if (handle) return 1;
  if (new_natoms != natoms) {
    fprintf(stderr, "number of atoms changed: %d -> %d\n", natoms, new_natoms);
    return 1;
  }
  close_file_read(handle);
  return 0;
}

#endif

///////////////////////////////////////////////////////////////////
//
// Plugin Interface
//
// ////////////////////////////////////////////////////////////////

static molfile_plugin_t desmond;
static molfile_plugin_t desmond_vel;

VMDPLUGIN_EXTERN int VMDPLUGIN_init (void) {
  /* Plugin for desmond trajectory files */
  ::memset(&desmond,0,sizeof(desmond));
  desmond.abiversion = vmdplugin_ABIVERSION;
  desmond.type = MOLFILE_PLUGIN_TYPE;
  desmond.name = "dtr";
  desmond.prettyname = "Desmond Trajectory";
  desmond.author = "D.E. Shaw Research";
  desmond.majorv = 2;
  desmond.minorv = 1;
  desmond.is_reentrant = VMDPLUGIN_THREADUNSAFE;

  desmond.filename_extension = "dtr,dtr/";
  desmond.open_file_read = open_file_read;
  desmond.read_timestep_metadata = read_timestep_metadata;
  desmond.read_next_timestep = read_next_timestep;
  desmond.close_file_read = close_file_read;

  desmond.open_file_write = open_file_write;
  desmond.write_structure = NULL;
  desmond.write_timestep = write_timestep;
  desmond.close_file_write = close_file_write;

  desmond_vel = desmond;
  desmond_vel.name = "dtrv";
  desmond_vel.prettyname = "Desmond Trajectory (with velocity)";

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  cb(v,(vmdplugin_t*)(void*)&desmond);
  cb(v,(vmdplugin_t*)(void*)&desmond_vel);
  //  cb(v,reinterpret_cast<vmdplugin_t*>(plugin_desmond));
  //  cb(v,reinterpret_cast<vmdplugin_t*>(plugin_desmond_vel));
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_EXTERN int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}


