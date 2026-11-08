"""Entry point for the frozen (PyInstaller) build. Same as the `pymol` script."""
import sys

import pymol

sys.exit(pymol.launch())
