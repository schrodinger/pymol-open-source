"""
Windows exceptions that Python handles must not abort the test run.

faulthandler's Windows handler also reports exceptions that are handled
afterwards (a DLL that fails to load and becomes an ImportError, or an
access violation that ctypes turns into OSError). Inside PyMOL on Python
3.15, dumping that report crashed the whole pytest session.
"""
import sys

import pytest

pytestmark = pytest.mark.skipif(sys.platform != 'win32',
                                reason='Windows structured exceptions')


def test_handled_access_violation():
    import ctypes

    # Same check as CPython's ctypes test_SEH: an invalid pointer
    # raises OSError instead of crashing the process.
    with pytest.raises(OSError, match='access violation'):
        ctypes.windll.kernel32.GetModuleHandleA(32)
