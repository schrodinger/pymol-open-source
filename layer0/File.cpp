/*
 * Copyright (c) Schrodinger, LLC.
 *
 * Basic file IO.
 */

#ifdef _WIN32
#include <vector>
#include <Windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "File.h"
#include "FileStream.h"
#include "MemoryDebug.h"

#include <fstream>
#include <string>

#include "pymol/zstring_view.h"

/**
 * Get the size from the current file pointer to the end of the file
 */
static long fgetsize(FILE *fp) {
  long filesize, current = ftell(fp);
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  fseek(fp, current, SEEK_SET);
  return filesize;
}

/**
 * Allocate memory and read the entire file from the given file pointer.
 * The file size is stored into the size pointer if not NULL.
 */
static char * fgetcontents(FILE *fp, long *size) {
  long filesize = fgetsize(fp);

  char *contents = pymol::malloc<char>(filesize + 255);
  if (!contents)
    return nullptr;

  if (1 != fread(contents, filesize, 1, fp)) {
    mfree(contents);
    return nullptr;
  }

  if (size)
    *size = filesize;

  contents[filesize] = '\0';
  return contents;
}

#ifdef _WIN32
FILE * pymol_fopen(const char * filename, const char * mode) {
  FILE *fp = fopen(filename, mode);

  if (!fp) {
    size_t len_filename = strlen(filename);
    std::vector<wchar_t> wfilename(len_filename + 1);
    std::vector<wchar_t> wmode(mode, mode + strlen(mode) + 1);

    if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS,
          filename, len_filename, wfilename.data(), wfilename.size()))
      return NULL;

    fp = _wfopen(wfilename.data(), wmode.data());
  }

  return fp;
}
#endif

/**
 * Allocate memory and read the entire file for the given filename.
 * The file size is stored into the size pointer if not NULL.
 */
char * FileGetContents(const char *filename, long *size) {
  char *contents;
  FILE *fp = pymol_fopen(filename, "rb");

  if (!fp)
    return nullptr;

  contents = fgetcontents(fp, size);
  fclose(fp);
  return contents;
}

namespace pymol
{
#ifdef _WIN32
std::wstring utf8_to_utf16(pymol::zstring_view utf8)
{
  std::wstring utf16(utf8.size(), L'\0');

  auto n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, utf8.c_str(), -1,
      &utf16[0], utf16.size() + 1);

  if (n == 0) {
    throw std::runtime_error("Could not convert to UTF-16");
  }

  utf16.resize(n - 1);
  return utf16;
}
#endif

static long istream_get_size(std::istream& file)
{
  auto current = file.tellg();
  file.seekg(0, std::ios::end);
  auto filesize = file.tellg();
  file.seekg(current);
  return filesize;
}

static std::string istream_get_contents(std::istream& file)
{
  auto filesize = istream_get_size(file);
  std::string contents(filesize, '\0');
  file.read(&contents[0], filesize);
  return contents;
}

std::string file_get_contents(pymol::zstring_view filename)
{
  std::ifstream file;
  fstream_open(file, filename, std::ios::in | std::ios::binary);
  return istream_get_contents(file);
}

} // namespace pymol
