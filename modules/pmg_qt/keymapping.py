'''
Module for translating Qt key codes to PyMOL key and "special" codes
'''

from pymol.Qt import QtCore
Qt = QtCore.Qt

DEBUG = False

keyMap = {
    Qt.Key_Escape: 27,
    Qt.Key_Tab: 9,
    Qt.Key_Backspace: 8,
    Qt.Key_Return: 13,
    Qt.Key_Enter: 13,
    Qt.Key_Delete: 127,
}

specialMap = {
    Qt.Key_Left: 100,
    Qt.Key_Up: 101,
    Qt.Key_Right: 102,
    Qt.Key_Down: 103,
    Qt.Key_PageUp: 104,
    Qt.Key_PageDown: 105,
    Qt.Key_Home: 106,
    Qt.Key_End: 107,
    Qt.Key_Insert: 108,
    Qt.Key_F1: 1,
    Qt.Key_F2: 2,
    Qt.Key_F3: 3,
    Qt.Key_F4: 4,
    Qt.Key_F5: 5,
    Qt.Key_F6: 6,
    Qt.Key_F7: 7,
    Qt.Key_F8: 8,
    Qt.Key_F9: 9,
    Qt.Key_F10: 10,
    Qt.Key_F11: 11,
    Qt.Key_F12: 12,
}


def get_modifiers(ev):
    '''Get modifers from event and translate into PyMOL modifier mask'''
    pymolmod = 0
    qtmodifiers = ev.modifiers()

    for mask, qtm in [
        (0x1, Qt.ShiftModifier),
        (0x2, Qt.MetaModifier),  # CTRL on Mac
        (0x2, Qt.ControlModifier),
        (0x4, Qt.AltModifier)
    ]:
        if qtmodifiers & qtm:
            pymolmod |= mask

    return pymolmod


def keyPressEventToPyMOLButtonArgs(ev):
    # translate modifier mask
    pymolmod = get_modifiers(ev)

    # Qt::Key_*
    key = ev.key()

    if key in specialMap:
        k = specialMap[key]

        # PyMOL_Special
        state = -2
    else:
        # PyMOL_Key
        state = -1

        # doesn't work on Mac: k = keyMap.get(key, ev.nativeVirtualKey())
        k = keyMap.get(key, -1)
        if k == -1:
            text = ev.text()
            if text:
                k = ord(text)

        # CTRL-<key>
        if k == -1 and (pymolmod & 0x2):
            k = key - 64

        # ALT-<key>
        if k != -1 and (pymolmod & 0x4):
            k = key

        if k > 255 or k < 0:
            if DEBUG:
                print('DEBUG: skipped: 0x%x 0x%x' % (key, k))
            return

    return (k, state, 0, 0, pymolmod)


def get_wheel_delta(ev):
    '''
    Get mouse wheel delta from event.
    Ignores horizontal scrolling (returns zero).
    '''
    try:
        # Qt4
        return ev.delta()
    except AttributeError:
        pass

    # Qt5
    angledelta = ev.angleDelta()
    delta_x = angledelta.x()
    delta_y = angledelta.y()

    if abs(delta_y) < abs(delta_x):
        # Shift+Wheel emulates horizontal scrolling
        if not (ev.modifiers() & Qt.ShiftModifier):
            return 0
        return delta_x

    return delta_y


def get_wheel_button(ev):
    '''
    Get mouse wheel button index (3 or 4) from event, or 0 if no vertial
    scrolling was detected.
    '''
    delta = get_wheel_delta(ev)
    if delta > 0:
        return 3
    if delta < 0:
        return 4
    return 0
