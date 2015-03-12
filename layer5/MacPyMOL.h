/*
 * MacPyMOL functions called from pymol repository
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifndef _H_MACPYMOL
#define _H_MACPYMOL

#ifdef _MACPYMOL_XCODE
void MacPyMOL_SetProgress(float value);
int MacPyMOL_doWindow(int code, int x,int y,int w, int h);
int MacPyMOL_fullScreenActive();
#endif

#ifdef _PYMOL_IOS
float IOS_getContentScaleFactor();
#define getContentScaleFactor IOS_getContentScaleFactor
#endif

#endif
