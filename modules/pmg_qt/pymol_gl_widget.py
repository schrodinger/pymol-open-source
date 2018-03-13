from __future__ import absolute_import

import sys
from pymol2 import SingletonPyMOL as PyMOL

from pymol.Qt import QtCore, QtOpenGL
from pymol.Qt import QtWidgets
Gesture = QtCore.QEvent.Gesture
Qt = QtCore.Qt

from .keymapping import get_modifiers

# don't import the heavy OpenGL (PyOpenGL) module
from pymol._cmd import glViewport


class PyMOLGLWidget(QtOpenGL.QGLWidget):
    '''
    PyMOL OpenGL Widget
    '''

    # mouse button map
    _buttonMap = {
        Qt.LeftButton: 0,
        Qt.MidButton: 1,
        Qt.RightButton: 2,
    }

    def __init__(self, parent):
        self.gui = parent

        # OpenGL context setup
        f = QtOpenGL.QGLFormat()
        f.setRgba(True)
        f.setDepth(True)
        f.setDoubleBuffer(True)

        from pymol.invocation import options

        # logic equivalent to layer5/main.cpp:launch

        if options.multisample:
            f.setSampleBuffers(True)

        if options.force_stereo != -1:
            # See layer1/Setting.h for stereo modes

            if options.stereo_mode in (0, 1, 12):
                # this effectively disables stereo detection
                # on Linux that is faulty in QGLWidget / PyQt5
                if not (options.stereo_mode == 0 and
                        sys.platform.startswith("linux")):
                    f.setStereo(True)

            if options.stereo_mode in (11, 12):
                f.setAccum(True)

            if options.stereo_mode in (0, 6, 7, 8, 9):
                f.setStencil(True)

        QtOpenGL.QGLWidget.__init__(self, f, parent=parent)

        if not self.isValid():
            raise RuntimeError('OpenGL initialization failed')

        f_actual = self.format()

        # report if quad buffer available
        options.stereo_capable = int(f_actual.stereo() or
                                     (options.force_stereo == 1))

        # feedback if stereo request failed
        if options.stereo_mode and (
                # QTBUG-59636 f.stereo() and not f_actual.stereo() or
                f.accum() and not f_actual.accum() or
                f.stencil() and not f_actual.stencil()):
            # cPyMOLGlobals_LaunchStatus_StereoFailed
            options.launch_status |= 0x1

        # feedback if multisample request failed
        if options.multisample and not f_actual.sampleBuffers():
            # cPyMOLGlobals_LaunchStatus_MultisampleFailed
            options.launch_status |= 0x2

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

    def mouseMoveEvent(self, ev):
        self.pymol.drag(int(self.fb_scale * ev.x()),
                        int(self.fb_scale * (self.height() - ev.y())),
                        get_modifiers(ev))

    def mousePressEvent(self, ev, state=0):
        if ev.button() not in self._buttonMap:
            return
        self.pymol.button(self._buttonMap[ev.button()], state,
                          int(self.fb_scale * ev.x()),
                          int(self.fb_scale * (self.height() - ev.y())),
                          get_modifiers(ev))

    def mouseReleaseEvent(self, ev):
        self.mousePressEvent(ev, 1)

    def wheelEvent(self, ev):
        pymolmod = get_modifiers(ev)
        try:
            delta = ev.delta()
        except AttributeError:
            # Qt5
            angledelta = ev.angleDelta()
            delta = angledelta.y()
            if abs(delta) < abs(angledelta.x()):
                # Shift+Wheel emulates horizontal scrolling
                if not (ev.modifiers() & Qt.ShiftModifier):
                    return
                delta = angledelta.x()
        if not delta:
            return
        button = 3 if delta > 0 else 4
        args = (int(self.fb_scale * ev.x()),
                int(self.fb_scale * (self.height() - ev.y())),
                pymolmod)
        self.pymol.button(button, 0, *args)
        self.pymol.button(button, 1, *args)

    ##########################
    # OpenGL
    ##########################

    def paintGL(self):
        glViewport(0, 0, int(self.fb_scale * self.width()),
                         int(self.fb_scale * self.height()))
        self.pymol.draw()
        self._timer.start(0)

    def resizeGL(self, w, h):
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
            self.updateFbScale(window)
            window.screenChanged.connect(self.updateFbScale)
        except AttributeError:
            # Fallback for Qt4
            self.fb_scale = 1.0

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
