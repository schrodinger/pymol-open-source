from pymol.Qt import QtGui


def textformat(color, style=''):
    fmt = QtGui.QTextCharFormat()
    fmt.setForeground(QtGui.QColor(color))

    for word in style.split():
        if word == 'bold':
            fmt.setFontWeight(QtGui.QFont.Bold)
        elif word == 'italic':
            fmt.setFontItalic(True)
        elif word.startswith('bg:'):
            fmt.setBackground(QtGui.QColor(word[3:]))
        else:
            print('unhandled style:', word)

    return fmt
