/*
OpenVR for PyMOL Copyright Notice
=====================================

The OpenVR for PyMOL source code is copyrighted, but you can freely use and
copy it as long as you don't change or remove any of the Copyright notices.
OpenVR for PyMOL is made available under the following open-source license
terms:

------------------------------------------------------------------------------
Copyright (c) 2018 EPAM Systems, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------

*/

#ifndef _H_OpenVRActionList
#define _H_OpenVRActionList

// system headers
#include <vector>

// local headers
#include "OpenVRAction.h"

struct OpenVRActionList {
  typedef std::vector<OpenVRAction*> ListType;
  ListType all;

  vr::VRActionSetHandle_t DefaultSet;

  OpenVRAction* LeftHand;
  OpenVRAction* RightHand;
  OpenVRAction* LGrip;
  OpenVRAction* RGrip;
  OpenVRAction* ToggleMenu;
  OpenVRAction* Laser;
  OpenVRAction* ActionSetNext;
  OpenVRAction* ActionSetPrev;
  OpenVRAction* Action1;
  OpenVRAction* Action2;
  OpenVRAction* Action3;

  explicit OpenVRActionList(vr::IVRInput* Input) {
    Input->GetActionSetHandle("/actions/pymol", &DefaultSet);

#define OPENVR_ADD_ACTION(set, Action, ...) all.push_back(Action = new OpenVRAction(Input, "/actions/" set "/in/" #Action, ##__VA_ARGS__));
    OPENVR_ADD_ACTION("pymol", LeftHand, OpenVRAction::TYPE_POSE);
    OPENVR_ADD_ACTION("pymol", RightHand, OpenVRAction::TYPE_POSE);
    OPENVR_ADD_ACTION("pymol", LGrip);
    OPENVR_ADD_ACTION("pymol", RGrip);
    OPENVR_ADD_ACTION("pymol", ToggleMenu);
    OPENVR_ADD_ACTION("pymol", Laser);
    OPENVR_ADD_ACTION("pymol", ActionSetNext);
    OPENVR_ADD_ACTION("pymol", ActionSetPrev);
    OPENVR_ADD_ACTION("pymol", Action1);
    OPENVR_ADD_ACTION("pymol", Action2);
    OPENVR_ADD_ACTION("pymol", Action3);
#undef OPENVR_ADD_ACTION
  }

  void Update(vr::IVRInput* Input) {
    vr::VRActiveActionSet_t activeSet = {DefaultSet};
    Input->UpdateActionState(&activeSet, sizeof(activeSet), 1);

    for (ListType::iterator i = all.begin(), iend = all.end(); i != iend; ++i) {
      (*i)->Update(Input);
    }
  }
};

#endif // _H_OpenVRActionList
