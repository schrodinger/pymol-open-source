import os
import sys
from pymol2 import SingletonPyMOL as PyMOL

import pymol

from pymol.Qt import QtCore
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets
Gesture = QtCore.QEvent.Gesture
Qt = QtCore.Qt

from .keymapping import get_modifiers
from .keymapping import get_wheel_button

# don't import the heavy OpenGL (PyOpenGL) module
from pymol._cmd import glViewport

# QOpenGLWidget is supposed to supersede QGLWidget, but has issues (e.g.
# no stereo support)
USE_QOPENGLWIDGET = int(
    os.getenv("PYMOL_USE_QOPENGLWIDGET") or
    (pymol.IS_MACOS and QtCore.QT_VERSION >= 0x50400))

if USE_QOPENGLWIDGET:
    BaseGLWidget = QtWidgets.QOpenGLWidget
    AUTO_DETECT_STEREO = False
else:
    from pymol.Qt import QtOpenGL
    BaseGLWidget = QtOpenGL.QGLWidget
    # only attempt stereo detection in Qt <= 5.6 (with 5.9+ on Linux I
    # get GL_DOUBLEBUFFER=0 with flickering when requesting stereo)
    AUTO_DETECT_STEREO = pymol.IS_WINDOWS or QtCore.QT_VERSION < 0x50700


class PyMOLGLWidget(BaseGLWidget):
    '''
    PyMOL OpenGL Widget

    Can be used as a context manager to make OpenGL context current.
    '''

    # mouse button map
    _buttonMap = {
        Qt.LeftButton: 0,
        Qt.MidButton: 1,
        Qt.RightButton: 2,
    }

    def __enter__(self):
        '''
        Context manager to make the OpenGL context "current".

        Fixes depth-buffer issue with QOpenGLWidget
        https://github.com/schrodinger/pymol-open-source/issues/25
        '''
        self.makeCurrent()
        pymol._cmd._pushValidContext(self.cmd._COb)

    def __exit__(self, exc_type, exc_value, traceback):
        pymol._cmd._popValidContext(self.cmd._COb)

    def __init__(self, parent):
        self.gui = parent
        self.fb_scale = 1.0

        # OpenGL context setup
        if USE_QOPENGLWIDGET:
            f = QtGui.QSurfaceFormat()
        else:
            f = QtOpenGL.QGLFormat()

        from pymol.invocation import options

        # logic equivalent to layer5/main.cpp:launch

        if options.multisample:
            f.setSamples(4)

        if options.force_stereo != -1:
            # See layer1/Setting.h for stereo modes

            if options.stereo_mode in (1, 12) or (
                    options.stereo_mode == 0 and AUTO_DETECT_STEREO):
                f.setStereo(True)

            if options.stereo_mode in (11, 12) and not USE_QOPENGLWIDGET:
                f.setAccum(True)

        if USE_QOPENGLWIDGET:
            super(PyMOLGLWidget, self).__init__(parent=parent)
            self.setFormat(f)
            self.setUpdateBehavior(QtWidgets.QOpenGLWidget.PartialUpdate)
        else:
            super(PyMOLGLWidget, self).__init__(f, parent=parent)

        # pymol instance
        self.pymol = PyMOL()
        self.pymol.start()
        self.cmd = self.pymol.cmd

        # capture python output for feedback
        import pcatch
        pcatch._install()

        # for passive move drag
        self.setMouseTracking(True)

        # for accepting keyboard input (command line, shortcuts)
        self.setFocusPolicy(Qt.ClickFocus)

        # for idle rendering
        self._timer = QtCore.QTimer()
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._pymolProcess)

        # drag n drop
        self.setAcceptDrops(True)

        # pinch-zoom
        self.grabGesture(Qt.PinchGesture)

    def sizeHint(self):
        # default 640 + internal_gui, 480 + internal_feedback
        return QtCore.QSize(860, 498)

    ##########################
    # Input Events
    ##########################

    def event(self, ev):
        if ev.type() == Gesture:
            return self.gestureEvent(ev)

        return super(PyMOLGLWidget, self).event(ev)

    def gestureEvent(self, ev):
        gesture = ev.gesture(Qt.PinchGesture)

        if gesture is None:
            return False

        if gesture.state() == Qt.GestureStarted:
            self.pinch_start_z = self.cmd.get_view()[11]

        changeFlags = gesture.changeFlags()

        if changeFlags & QtWidgets.QPinchGesture.RotationAngleChanged:
            delta = gesture.lastRotationAngle() - gesture.rotationAngle()
            self.cmd.turn('z', delta)

        if changeFlags & QtWidgets.QPinchGesture.ScaleFactorChanged:
            view = list(self.cmd.get_view())

            # best guess for https://bugreports.qt.io/browse/QTBUG-48138
            totalscalefactor = gesture.totalScaleFactor()
            if totalscalefactor == 1.0:
                totalscalefactor = gesture.scaleFactor()

            z = self.pinch_start_z / totalscalefactor
            delta = z - view[11]
            view[11] = z
            view[15] -= delta
            view[16] -= delta
            self.cmd.set_view(view)

        return True

    def _event_x_y_mod(self, ev):
        return (
            int(self.fb_scale * ev.x()),
            int(self.fb_scale * (self.height() - ev.y())),
            get_modifiers(ev),
        )

    def mouseMoveEvent(self, ev):
        self.pymol.drag(*self._event_x_y_mod(ev))

    def mousePressEvent(self, ev, state=0):
        if ev.button() not in self._buttonMap:
            return
        self.pymol.button(self._buttonMap[ev.button()], state,
                          *self._event_x_y_mod(ev))

    def mouseReleaseEvent(self, ev):
        self.mousePressEvent(ev, 1)

    def wheelEvent(self, ev):
        button = get_wheel_button(ev)
        if not button:
            return
        args = self._event_x_y_mod(ev)
        self.pymol.button(button, 0, *args)
        self.pymol.button(button, 1, *args)

    ##########################
    # OpenGL
    ##########################

    def paintGL(self):
        if not USE_QOPENGLWIDGET:
            glViewport(0, 0, int(self.fb_scale * self.width()),
                         int(self.fb_scale * self.height()))
        self.pymol.draw()
        self._timer.start(0)

    def resizeGL(self, w, h):
        if USE_QOPENGLWIDGET:
            w = int(w * self.fb_scale)
            h = int(h * self.fb_scale)

        self.pymol.reshape(w, h, True)

    def updateFbScale(self, context):
        '''Update PyMOL's display scale factor from the window or screen context
        @type context: QWindow or QScreen
        '''
        self.fb_scale = context.devicePixelRatio()
        try:
            self.cmd.set('display_scale_factor', int(self.fb_scale))
        except BaseException as e:
            # fails with modal draw (mpng ..., modal=1)
            print(e)

    def initializeGL(self):
        # Scale framebuffer for Retina displays
        try:
            window = self.windowHandle()

            # QOpenGLWidget workaround
            if window is None:
                window = self.parent().windowHandle()

            self.updateFbScale(window)
            window.screenChanged.connect(self.updateFbScale)
            window.screen().physicalDotsPerInchChanged.connect(
                    lambda dpi: self.updateFbScale(window))

        except AttributeError:
            # Fallback for Qt4
            pass

    def _pymolProcess(self):
        idle = self.pymol.idle()
        if idle or self.pymol.getRedisplay():
            self.update()

        self._timer.start(20)

    ##########################
    # drag n drop
    ##########################

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    url = url.toLocalFile()
                else:
                    url = url.toString()
                self.gui.load_dialog(url)
            event.accept()
