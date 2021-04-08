from pymol.Qt import *


class UpdateLock:
    """
    Locking mechanism to prevent circular signal/slot updates.

    Decorator notation:

        >>> @lock.skipIfCircular
        >>> def updatesomething():
        ...     dosomething()

    Context manager notation:

        >>> def updatesomething():
        ...     with lock:
        ...         lock.acquire()  # exits block if already acquired
        ...         dosomething()

    """

    class LockFailed(Exception):
        pass

    def __init__(self, silent_exc_types=()):
        self.primed = False
        self.acquired = False
        self.silent_exc_types = tuple(silent_exc_types)

    def __enter__(self):
        assert not self.primed, 'missing acquire()'
        self.primed = True

    def __exit__(self, exc_type, exc_val, exc_tb):
        assert not self.primed, 'missing acquire()'

        if exc_type == self.LockFailed:
            return True

        assert self.acquired, 'inconsistency!?'
        self.acquired = False

        if exc_type in self.silent_exc_types:
            return True

    def acquire(self):
        assert self.primed, 'missing with ...():'
        self.primed = False

        if self.acquired:
            raise self.LockFailed

        self.acquired = True

    def skipIfCircular(self, func):
        def wrapper(*args, **kwargs):
            with self:
                self.acquire()
                return func(*args, **kwargs)
        return wrapper


class WidgetMenu(QtWidgets.QMenu):
    """
    QMenu that represents a single widget that pops up under a push button.

    >>> btn = QtWidgets.QPushButton()
    >>> btn.setMenu(WidgetMenu().setSetupUi(setupUi))
    """

    def focusNextPrevChild(self, next):
        '''Overload which prevents menu-like tab action'''
        return QtWidgets.QWidget.focusNextPrevChild(self, next)

    def setWidget(self, widget):
        self.clear()
        action = QtWidgets.QWidgetAction(self)
        action.setDefaultWidget(widget)
        self.addAction(action)
        return self

    def setSetupUi(self, setupUi):
        '''Use a setup function for the widget. The widget will be created
        and initialized as "setupUi(widget)" before the menu is shown for
        the first time.'''

        @self.aboutToShow.connect
        def _():
            self.aboutToShow.disconnect()
            widget = QtWidgets.QWidget()
            setupUi(widget)
            self.setWidget(widget)

        return self


class AsyncFunc(QtCore.QThread):
    """
    Decorator to call a function asynchronous.

    Signals:
    - returned(result)
        Emitted if the function successfully returned
    - finished(result, exception)
        Emitted on success or failure. On success, exception is None,
        on failure result is None.

    Example:
    >>> asyncsquare = AsyncFunc(lambda x: x * x, print)
    >>> asyncsquare(5)
    25

    """
    # Warning: PySide crashes if passing None to an "object" type signal
    returned = QtCore.Signal(object)
    finished = QtCore.Signal(tuple)

    def __init__(self, func, returnslot=None, finishslot=None):
        super(AsyncFunc, self).__init__()
        self.func = func
        if returnslot is not None:
            self.returned.connect(returnslot)
        if finishslot is not None:
            self.finished.connect(finishslot)

    def __call__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.start()

    def run(self):
        result = None
        exception = None
        try:
            result = self.func(*self.args, **self.kwargs)
            self.returned.emit(result)
        except Exception as e:
            exception = e
        except:
            exception = Exception()
        self.finished.emit((result, exception))


class MainThreadCaller(QtCore.QObject):
    """
    Allows calling a GUI function from a non-main thread. Will pause
    the current thread until the function has returned from the main
    thread.

        >>> # in main thread
        >>> callInMainThread = MainThreadCaller()

        >>> # in async thread
        >>> callInMainThread(lambda: 123)
        123

    Note: QMetaObject.invokeMethod with BlockingQueuedConnection could
    potentially be used to achieve the same goal.
    """
    mainthreadrequested = QtCore.Signal(object)

    RESULT_RETURN = 0
    RESULT_EXCEPTION = 1

    def __init__(self):
        super(MainThreadCaller, self).__init__()
        self.waitcondition = QtCore.QWaitCondition()
        self.mutex = QtCore.QMutex()
        self.results = {}
        self.mainthreadrequested.connect(self._mainThreadAction)

    def _mainThreadAction(self, func):
        try:
            self.results[func] = (self.RESULT_RETURN, func())
        except Exception as ex:
            self.results[func] = (self.RESULT_EXCEPTION, ex)
        self.waitcondition.wakeAll()

    def __call__(self, func):
        if self.thread() is QtCore.QThread.currentThread():
            return func()

        self.mainthreadrequested.emit(func)

        while True:
            self.mutex.lock()
            self.waitcondition.wait(self.mutex)
            self.mutex.unlock()
            try:
                result = self.results.pop(func)
                break
            except KeyError:
                print(type(self).__name__ + ': result was not ready')

        if result[0] == self.RESULT_EXCEPTION:
            raise result[1]

        return result[1]


def connectFontContextMenu(widget):
    """
    Connects a custom context menu with a "Select Font..." entry
    to the given widget.

    @type widget: QWidget
    """
    widget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)

    @widget.customContextMenuRequested.connect
    def _(pt):
        menu = widget.createStandardContextMenu()
        menu.addSeparator()
        action = menu.addAction("Select Font...")

        @action.triggered.connect
        def _():
            font, ok = QtWidgets.QFontDialog.getFont(widget.font(), widget,
                    "Select Font", QtWidgets.QFontDialog.DontUseNativeDialog)
            if ok:
                widget.setFont(font)

        menu.exec_(widget.mapToGlobal(pt))


def getSaveFileNameWithExt(*args, **kwargs):
    """
    Return a file name, append extension from filter if no extension provided.
    """
    import os, re

    fname, filter = QtWidgets.QFileDialog.getSaveFileName(*args, **kwargs)

    if not fname:
        return ''

    if '.' not in os.path.split(fname)[-1]:
        m = re.search(r'\*(\.[\w\.]+)', filter)
        if m:
            # append first extension from filter
            fname += m.group(1)

    return fname


def getMonospaceFont(size=9):
    """
    Get the best looking monospace font for the current platform
    """
    import sys

    if sys.platform == 'darwin':
        family = 'Monaco'
        size += 3
    elif sys.platform == 'win32':
        family = 'Consolas'
    else:
        family = 'Monospace'

    font = QtGui.QFont(family, size)
    font.setStyleHint(font.Monospace)

    return font


def loadUi(uifile, widget):
    """
    Load .ui file into widget

    @param uifile: filename
    @type uifile: str
    @type widget: QtWidgets.QWidget
    """
    if PYQT_NAME.startswith('PyQt'):
        m = __import__(PYQT_NAME + '.uic')
        return m.uic.loadUi(uifile, widget)
    elif PYQT_NAME == 'PySide2':
        try:
            import pyside2uic as pysideuic
        except ImportError:
            pysideuic = None
    else:
        import pysideuic

    if pysideuic is None:
        import subprocess
        p = subprocess.Popen(['uic', '-g', 'python', uifile],
                             stdout=subprocess.PIPE)
        source = p.communicate()[0]
        # workaround for empty retranslateUi bug
        source += b'\n' + b' ' * 8 + b'pass'
    else:
        import io
        stream = io.StringIO()
        pysideuic.compileUi(uifile, stream)
        source = stream.getvalue()

    ns_locals = {}
    exec(source, ns_locals)

    if 'Ui_Form' in ns_locals:
        form = ns_locals['Ui_Form']()
    else:
        form = ns_locals['Ui_Dialog']()

    form.setupUi(widget)
    return form


class PopupOnException:
    """
    Context manager which shows a message box if an exception is raised.

        >>> with PopupOnException():
        ...     # do something

    Decorator support:

        >>> @PopupOnException.decorator
        ... def foo():
        ...     # do something

    """

    @classmethod
    def decorator(cls, func):
        def wrapper(*args, **kwargs):
            with cls():
                return func(*args, **kwargs)
        return wrapper

    def __enter__(self):
        pass

    def __exit__(self, exc_type, e, tb):
        if e is not None:
            import traceback
            QMB = QtWidgets.QMessageBox

            parent = QtWidgets.QApplication.focusWidget()

            msg = str(e) or 'unknown error'
            msgbox = QMB(QMB.Critical, 'Error', msg, QMB.Close, parent)
            msgbox.setDetailedText(''.join(traceback.format_tb(tb)))
            msgbox.exec_()

        return True


def conda_ask_install(packagespec, channel=None, msg="", parent=None, url=""):
    """
    Install a conda package
    """
    return True
