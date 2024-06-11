#!/usr/bin/env python3
import re


print(re.findall(r'_PyMOL_VERSION "(.*)"', open('layer0/Version.h').read())[0])
