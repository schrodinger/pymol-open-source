# NOTE: This only works for module-based PyMOL builds.
# Embedded versions of PyMOL call PyUnicode_SetDefaultEncoding at startup
import sys
sys.setdefaultencoding("utf-8")
