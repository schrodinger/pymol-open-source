/*
 * MacPyMOL functions called from pymol repository
 *
 * (c) 2014 Schrodinger, Inc.
 */

#ifdef _MACPYMOL_XCODE
#ifndef _H_MACPYMOL
#define _H_MACPYMOL

void MacPyMOLSetProgress(float value);
int MacPyMOL_doWindow(int code, int x,int y,int w, int h);
int MacPyMOL_fullScreenActive();

#endif
#endif
