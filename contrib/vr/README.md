# OpenVR extension for PyMOL

Copyright (c) 2018 [EPAM Systems, Inc.](https://www.epam.com/)

Published under an MIT licence, see [LICENSE](LICENSE).

Original release URL:
https://github.com/epam/pymol-open-source/tree/openvr


PyMOL OpenVR Requirements
=========================

1. In order to experience VR in PyMOL you should have a VR system installed and ready. Currently we support HTC Vive only. Please visit VIVE Setup site (https://www.vive.com/us/setup/vive/) and follow instructions to install software and setup your room.

2. You should have at least SteamVR 1533664367 installed (August 7, 2018), however, it is recommended to have the most recent Steam version (installed during the Vive setup). Different versions of OpenVR/SteamVR might have different API (a set of features), so if you have an outdated Steam the chance is that the application will not run.

3. Test your PC (https://dl.vive.com/oobe/ViveCheck.exe) to be sure the performance is good enough. Recommended computer spec listed on HTC Vive website are:

- Processor: Intel™ Core™ i5-4590 or AMD FX™ 8350, equivalent or better.
- Graphics: NVIDIA GeForce™ GTX 1060 or AMD Radeon™ RX 480, equivalent or better.
- Memory: 4 GB RAM or more.
- Video output: 1x HDMI 1.4 port, or DisplayPort 1.2 or newer.
- USB: 1x USB 2.0 port or newer.
- Operating system: Windows™ 7 SP1, Windows™ 8.1 or later or Windows™ 10.

Note: You can still install and run PyMOL on a system without VR hardware. In this case you cannot use OpenVR stereo mode but the rest of PyMOL should function perfectly.
