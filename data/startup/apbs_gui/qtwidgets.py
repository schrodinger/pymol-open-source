from pymol.Qt import QtCore, QtWidgets


class ResizableMessageBox(QtWidgets.QMessageBox):

    _EVENT_TYPES = (
            QtCore.QEvent.UpdateRequest,
            QtCore.QEvent.WinIdChange,
            QtCore.QEvent.ShowToParent,
            )

    _UNWANTED_WINDOW_FLAGS = (
            QtCore.Qt.MSWindowsFixedSizeDialogHint |
            0)

    def _make_resizable(self):
        textEdit = self.findChild(QtWidgets.QTextEdit)
        if textEdit is None:
            return

        self.setSizeGripEnabled(True)

        ex = QtWidgets.QSizePolicy.Expanding
        for w in [self, textEdit]:
            w.setMaximumSize(0xffffff, 0xffffff)
            w.setSizePolicy(ex, ex)

    def event(self, e):
        if e.type() in self._EVENT_TYPES:
            self._make_resizable()

        return super(ResizableMessageBox, self).event(e)
