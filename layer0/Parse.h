

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
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_Parse
#define _H_Parse

const char *ParseWordCopy(char *dst, const char *src, int n);
const char *ParseWordNumberCopy(char *dst, const char *src, int n);
const char *ParseWord(char *dst, const char *src, int n);
const char *ParseNTrim(char *q, const char *p, int n);
const char *ParseNTrimRight(char *q, const char *p, int n);
const char *ParseNSkip(const char *p, int n);
const char *ParseCommaCopy(char *q, const char *p, int n);
const char *ParseSkipEquals(const char *p);
const char *ParseIntCopy(char *q, const char *p, int n);
const char *ParseAlphaCopy(char *q, const char *p, int n);
int ParseFloat3List(const char *p, float *vals);
const char *ParseNextLine(const char *p);
const char *ParseNCopy(char *dst, const char *src, int n);

/*
 * non-const overloads
 */
inline char *ParseWordCopy(char *q, char *p, int n) {
  return const_cast<char*>(ParseWordCopy(q, const_cast<const char*>(p), n));
}
inline char *ParseWordNumberCopy(char *q, char *p, int n) {
  return const_cast<char*>(ParseWordNumberCopy(q, const_cast<const char*>(p), n));
}
inline char *ParseWord(char *q, char *p, int n) {
  return const_cast<char*>(ParseWord(q, const_cast<const char*>(p), n));
}
inline char *ParseNTrim(char *q, char *p, int n) {
  return const_cast<char*>(ParseNTrim(q, const_cast<const char*>(p), n));
}
inline char *ParseNTrimRight(char *q, char *p, int n) {
  return const_cast<char*>(ParseNTrimRight(q, const_cast<const char*>(p), n));
}
inline char *ParseNSkip(char *p, int n) {
  return const_cast<char*>(ParseNSkip(const_cast<const char*>(p), n));
}
inline char *ParseCommaCopy(char *q, char *p, int n) {
  return const_cast<char*>(ParseCommaCopy(q, const_cast<const char*>(p), n));
}
inline char *ParseSkipEquals(char *p) {
  return const_cast<char*>(ParseSkipEquals(const_cast<const char*>(p)));
}
inline char *ParseIntCopy(char *q, char *p, int n) {
  return const_cast<char*>(ParseIntCopy(q, const_cast<const char*>(p), n));
}
inline char *ParseAlphaCopy(char *q, char *p, int n) {
  return const_cast<char*>(ParseAlphaCopy(q, const_cast<const char*>(p), n));
}
inline char *ParseNextLine(char *p) {
  return const_cast<char*>(ParseNextLine(const_cast<const char*>(p)));
}
inline char *ParseNCopy(char *q, char *p, int n) {
  return const_cast<char*>(ParseNCopy(q, const_cast<const char*>(p), n));
}

#endif
