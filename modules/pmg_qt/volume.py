"""
Volume Color Map Editor Panel.
"""

import itertools, math

import pymol

from pymol.Qt import QtGui, QtCore
from pymol.Qt import QtWidgets

Qt = QtCore.Qt

DOT_RADIUS = 5
ALPHA_LOG_BASE = 10.0

DEFAULT_COLORS = [
    (1., 1., 0.),
    (1., 0., 0.),
    (0., 0., 1.),
    (0., 1., 0.),
    (0., 1., 1.),
    (1., 0., 1.),
]

EPS = 1e-6

DEFAULT_TEXT_DIALOG_WIDTH = 500

VOLUME_HELP = '''
VOLUME PANEL HELP

--------------------------------------------------
Canvas Mouse Actions (no Point under Cursor)

  L-Click            Add point
  CTRL+L-Click       Add 3 points (isosurface)

  CTRL+R-Drag        Zoom in

--------------------------------------------------
Mouse Actions with Point under Cursor

  L-Click            Edit point color
  R-Click            Edit point value
  SHIFT+R-Click      Edit point opacity
  CTRL+L-Click       Edit color of 3 points

  M-Click            Remove Point
  SHIFT+L-Click      Remove Point
  CTRL+M-Click       Remove 3 points
  CTRL+SHIFT+L-Click Remove 3 points

  L-Drag             Move point
  CTRL+L-Drag        Move 3 points (horizontal only)
  R-Drag             Move point along one axis only

--------------------------------------------------
L = Left mouse button
M = Middle mouse button
R = Right mouse button

--------------------------------------------------
See also the "volume_color" command for getting and
setting volume colors on the command line.
'''


class VolumeEditorWidget(QtWidgets.QWidget):
    def __init__(self, parent=None, volume_name='', cmd=None):
        super(VolumeEditorWidget, self).__init__(parent)
        self.setObjectName("volume_editor_widget")
        self.setMouseTracking(True)
        self.points = []
        self.point = -1
        self.color_cycle = itertools.cycle(DEFAULT_COLORS)
        self.path = None
        self.cmd = cmd
        self.volume_name = volume_name
        self.real_time = True
        self.dragged = False
        self.vmin = 0.0
        self.vmax = 1.0
        self.original_vmin = 0.0
        self.original_vmax = 1.0
        self.amax = 1.0
        self.left_margin = 35
        self.bottom_margin = 20
        self.constrain = None
        self.init_pos = None
        self.zoom_pos = None
        self.ignore_set_colors = False
        self.color_dialog = None
        self.text_dialog = None
        self.text_boxes = {}
        self.hover_point = -1

    def sizeHint(self):
        return QtCore.QSize(600, 200)

    def paintGrid(self, painter, rect):
        pen = QtGui.QPen(self.line_color)
        painter.setPen(pen)
        x0 = rect.x()
        x1 = rect.x() + rect.width()
        y0 = rect.y()
        y1 = rect.y() + rect.height()
        painter.drawLine(x0, y1, x1, y1)
        painter.drawLine(x0, y0, x0, y1)
        h = rect.height()
        num_lines = 10
        pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        for line in range(1, num_lines):
            y = y0 + h * (1.0 - self.alphaToY(line / float(num_lines)))
            painter.drawLine(x0, y, x1, y)

    def paintColorDots(self, painter, rect):
        pen = QtGui.QPen(Qt.gray)
        pen.setStyle(Qt.SolidLine)
        painter.setPen(pen)
        scaled_pts = []
        h = rect.height()
        for point in self.points:
            x, y, r, g, b = point
            x = rect.left() + rect.width() * self.dataToX(x)
            y = rect.top() + rect.height() * (1.0 - self.alphaToY(y))
            if scaled_pts:
                painter.drawLine(scaled_pts[-1][0], scaled_pts[-1][1], x, y)
            scaled_pts.append((x, y, r, g, b))

        for x, y, r, g, b in scaled_pts:
            painter.setBrush(QtGui.QColor(255 * r, 255 * g, 255 * b))
            painter.drawEllipse(x - DOT_RADIUS, y - DOT_RADIUS, 2 * DOT_RADIUS,
                                2 * DOT_RADIUS)

        if 0 <= self.hover_point < len(scaled_pts):
            # use larger radius for hover dot
            radius = DOT_RADIUS + 2
            x, y, r, g, b = scaled_pts[self.hover_point]
            painter.setBrush(QtGui.QColor(255 * r, 255 * g, 255 * b))
            painter.drawEllipse(x - radius, y - radius, 2 * radius,
                                2 * radius)

    def paintHistogram(self, painter, rect):
        if self.path:
            vrange = self.original_vmax - self.original_vmin
            if vrange == 0.0:
                return
            norm_min = (self.vmin - self.original_vmin) / vrange
            norm_max = (self.vmax - self.original_vmin) / vrange
            h = rect.height() - 2
            dnorm = norm_max-norm_min
            iwidth = 1.0 / (rect.width())
            painter_path = QtGui.QPainterPath()
            for i in range(rect.width()):
                pos = (i * iwidth * dnorm + norm_min) * len(self.path)
                ipos = int(pos)
                if pos < 0 or ipos >= len(self.path) - 1:
                    continue
                y0 = self.path[ipos][1]
                y1 = self.path[ipos+1][1]
                y = y0 + (y1-y0) * (pos - int(pos)) # lerp
                x = rect.left() + i
                y = h - self.alphaToY(y) * h + 1
                if painter_path.elementCount() == 0:
                    painter_path.moveTo(x, y)
                else:
                    painter_path.lineTo(x, y)

            pen = QtGui.QPen(Qt.red)
            pen.setStyle(Qt.SolidLine)
            painter.setPen(pen)
            painter.drawPath(painter_path)

    def paintValueBox(self,
                      painter,
                      font_metrics,
                      x,
                      y,
                      right_just,
                      value,
                      format="%.3f"):
        s = format % value
        sw = font_metrics.width(s)
        sh = font_metrics.height()
        if right_just:
            rect = QtCore.QRect(x - sw - 4, y - sh, sw + 4, sh + 2)
        else:
            rect = QtCore.QRect(x, y - sh, sw + 4, sh + 2)
        painter.fillRect(rect,
                QtGui.QColor(96, 96, 128) if self.line_color == Qt.lightGray else
                QtGui.QColor(0xFF, 0xFF, 0xFF))
        painter.drawRect(rect)
        painter.drawText(rect.x() + 2, y - 2, s)
        return rect

    def paintAxes(self, painter, rect):
        low = int(math.ceil(self.vmin))
        hi = int(math.floor(self.vmax)) + 1
        pen = painter.pen()
        pen.setStyle(Qt.SolidLine)
        pen.setColor(self.line_color)
        painter.setPen(pen)
        fm = QtGui.QFontMetrics(painter.font())
        fw = fm.averageCharWidth()
        fh = fm.height() - 2
        x0 = rect.left()
        x1 = rect.right()
        y0 = rect.bottom()
        y1 = y0 + 4
        lastx = x0

        # horizontal axis
        for tick in range(low, hi):
            s = str(tick)
            w = fw * len(s)
            x = x0 + w / 2 + rect.width() * (tick - self.vmin
                                             ) / float(self.vmax - self.vmin)
            if x - lastx > w + 2 * fw:
                painter.drawLine(x, y0, x, y1)
                painter.drawText(x - w / 2, y1 + fh - 2, s)
                lastx = x

        #vertical axis
        x1 = rect.left()
        lasty = y0
        for tick in range(1, 10):
            t = tick / 10.0
            y = y0 - self.alphaToY(t) * rect.height()
            if lasty - y > fh and y > 2 * fh:
                painter.drawLine(x1 - 5, y, x1, y)
                painter.drawText(x1 - 5 - 3 * fw, y - 2 + fh / 2, str(t))
                lasty = y

        # text boxes
        self.text_boxes["vmin"] = self.paintValueBox(painter, fm,
                                                     rect.left(), y1 + fh,
                                                     False, self.vmin)
        self.text_boxes["vmax"] = self.paintValueBox(painter, fm,
                                                     rect.right(), y1 + fh,
                                                     True, self.vmax)
        self.text_boxes["amax"] = self.paintValueBox(
            painter, fm, x0 - 4 * fw, 2 + fh, False, self.amax, format="%.2f")

    def paintZoomArea(self, painter, rect):
        if self.init_pos and self.zoom_pos:
            rect.setLeft(self.init_pos.x())
            rect.setRight(self.zoom_pos.x())
            painter.fillRect(rect, QtGui.QBrush(QtGui.QColor(0, 64, 128, 128)))

    def paintEvent(self, event):
        """
        Paints the editor widget.
        """
        painter = QtGui.QPainter()
        self.paint_rect = event.rect()
        self.paint_rect.adjust(self.left_margin, 0, 0, -self.bottom_margin)

        # tweak color depening on the panel floating state
        # disabled: always use default style
        is_floating = True  # self.parent().parent().isFloating()
        self.line_color = Qt.darkGray if is_floating else Qt.lightGray

        painter.begin(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        self.paintGrid(painter, self.paint_rect)
        self.paintAxes(painter, self.paint_rect)
        painter.setClipRect(self.paint_rect)
        self.paintHistogram(painter, self.paint_rect)
        painter.setClipping(False)
        self.paintColorDots(painter, self.paint_rect)
        self.paintZoomArea(painter, self.paint_rect)
        painter.end()

    def enterValue(self, title, value, min_value, max_value):
        """
        Handles entering new values into alpha / min / max boxes.
        """
        new_value, status = QtWidgets.QInputDialog.getDouble(
            self, "", title, value, min_value, max_value, decimals=6)
        if status:
            return new_value

        return value

    def mousePressEvent(self, event):
        # process textbox clicks
        if event.button() == Qt.LeftButton:
            for key, rect in self.text_boxes.items():
                if rect.contains(event.pos()):
                    if key == "amax":
                        self.amax = self.enterValue("Maximum Alpha Value",
                                                    self.amax, EPS, 1.0)
                    elif key == "vmin":
                        self.vmin = self.enterValue("Minimum Data Value",
                                                    self.vmin, -1e8,
                                                    self.vmax - EPS)
                    else:
                        self.vmax = self.enterValue("Maximum Data Value",
                                                    self.vmax, self.vmin + EPS,
                                                    1e8)
                    self.repaint()
                    return

        if self.paint_rect.adjusted(
            -DOT_RADIUS, -DOT_RADIUS, DOT_RADIUS, DOT_RADIUS).contains(
            event.pos()):
            self.dragged = False
            self.point = self.findPoint(event.pos())
            self.init_pos = event.pos()
            self.zoom_pos = None
            self.constraint = None
            if self.point < 0 and event.button() == Qt.LeftButton:
                self.addPoint(
                    event.pos(), event.modifiers() == Qt.ControlModifier)
                # suppress color picker
                self.dragged = True

    def mouseReleaseEvent(self, event):
        if not self.dragged and self.point >= 0:
            if event.button() == Qt.RightButton:
                x, y, r, g, b = self.points[self.point]
                # in 2.0: help says NoModifier, implemented is ControlModifier
                if event.modifiers() in (Qt.ControlModifier, Qt.NoModifier):
                    value = self.points[self.point][0]
                    prev_x = self.points[self.point-1][0] if self.point > 0 else self.vmin
                    next_x = self.points[self.point+1][0] if self.point < len(self.points)-1 else self.vmax
                    x = self.enterValue("Data value",
                        value, prev_x, next_x)
                elif event.modifiers() == Qt.ShiftModifier:
                    value = self.points[self.point][1]
                    y = self.enterValue("Alpha value (opacity)",
                        value, 0.0, 1.0)
                self.points[self.point] = (x, y, r, g, b)
                self.repaint()
                if self.real_time:
                    self.updateVolumeColors()
            if (event.button() == Qt.MidButton or
                (event.button() == Qt.LeftButton and
                 event.modifiers() & Qt.ShiftModifier)):
                self.removePoints(
                    event.modifiers() & Qt.ControlModifier)
            elif event.button() == Qt.LeftButton:
                self.setPointColor(self.point,
                        event.modifiers() == Qt.ControlModifier)

        self.point = -1
        self.hover_point = -1

        if self.init_pos and self.zoom_pos:
            if self.init_pos.x() != self.zoom_pos.x():
                # zoom in
                self.vmin, self.vmax = sorted([
                    self.xToData(self.convertX(self.init_pos.x())),
                    self.xToData(self.convertX(self.zoom_pos.x()))])

        self.zoom_pos = None
        self.repaint()
        self.updateVolumeColors()

    def changePointColor(self, color):
        """
        Changes color of specified point or three points.
        """
        x, y, _, _, _ = self.points[self.color_point]
        r = color.redF()
        g = color.greenF()
        b = color.blueF()
        self.points[self.color_point] = (x, y, r, g, b)
        if self.color_triple:
            # set color for all three points
            if self.color_point > 0:
                x, y, _, _, _ = self.points[self.color_point - 1]
                self.points[self.color_point - 1] = (x, y, r, g, b)
            if self.color_point < len(self.points) - 1:
                x, y, _, _, _ = self.points[self.color_point + 1]
                self.points[self.color_point + 1] = (x, y, r, g, b)

        self.repaint()

    def updatePointColor(self, color):
        """
        This is called when color is changed in real time.
        """
        self.changePointColor(color)
        if self.real_time:
            self.updateVolumeColors()

    def colorDialogClosed(self, result):
        """
        This is called when color dialog is closed.
        """
        if result == QtWidgets.QDialog.Accepted:
            color = self.color_dialog.currentColor()
        else:
            color = self.original_color
        self.changePointColor(color)
        self.updateVolumeColors()

    def setPointColor(self, point, triple):
        """
        Opens color picker and sets color of one or three points.
        """
        self.color_point = point
        self.color_triple = triple
        _, _, r, g, b = self.points[self.color_point]
        if not self.color_dialog:
            self.color_dialog = QtWidgets.QColorDialog(self)
            self.color_dialog.currentColorChanged.connect(
                self.updatePointColor)
            self.color_dialog.finished.connect(self.colorDialogClosed)
        self.original_color = QtGui.QColor(255 * r, 255 * g, 255 * b)
        self.color_dialog.setCurrentColor(self.original_color)
        # open modal color dialog
        self.color_dialog.open()

    def toggleRealTimeUpdates(self, value):
        self.real_time = value

    def getColors(self):
        """
        Returns color map colors and coordinates as a flat list.

        @rtype: list of float
        @return: List of x, r, g, b, y values.
        """
        colors = []
        for point in self.points:
            x, y, r, g, b = point
            colors.append(x)
            colors.append(r)
            colors.append(g)
            colors.append(b)
            colors.append(y)
        return colors

    def convertX(self, x_pos):
        """
        Converts mouse X position within the widget to <0,1> range.
        """
        x = (x_pos - self.left_margin) / float(self.width() - self.left_margin)
        return min(max(x, 0.0), 1.0)

    def convertY(self, y_pos):
        """
        Converts mouse Y position within the widget to <0,1> range.
        """
        y = 1.0 - y_pos / float(self.height() - self.bottom_margin)
        return min(max(y, 0.0), 1.0)

    def xToData(self, x):
        """
        Converts <0, 1> value to actual data range.
        """
        return self.vmin + x * (self.vmax - self.vmin)

    def dataToX(self, d):
        """
        Converts data value to <0, 1> range.
        """
        return (d - self.vmin) / (self.vmax - self.vmin)

    def yToAlpha(self, y):
        """
        Converts <0, 1> normalized Y position to alpha value.
        """
        y = (ALPHA_LOG_BASE**y - 1.0) / (ALPHA_LOG_BASE - 1.0)
        return y * self.amax

    def alphaToY(self, a):
        """
        Converts alpha value to <0, 1> normalized Y position.
        """
        if self.amax == 0.0:
            return 0.0
        y = a / self.amax
        return math.log(1.0 + (ALPHA_LOG_BASE - 1.0) * y, ALPHA_LOG_BASE)

    def updateVolumeColors(self):
        """
        Updates volume colors in PyMOL display.
        """
        self.ignore_set_colors = True
        self.cmd.volume_color(self.volume_name, self.getColors())
        self.ignore_set_colors = False

    def mouseMoveEvent(self, event):
        if event.buttons() in (Qt.LeftButton, Qt.RightButton):
            if (event.buttons() == Qt.RightButton and
                    event.modifiers() == Qt.ControlModifier):
                # zoom in
                self.zoom_pos = event.pos()
                self.repaint()
            elif self.point >= 0:
                self.dragged = True
                if event.buttons(
                ) & Qt.RightButton and not self.constraint:
                    # constrained movement
                    dpos = event.pos() - self.init_pos
                    self.constraint = 'x' if (
                        abs(dpos.x()) > abs(dpos.y())) else 'y'
                self.movePoints(event)
        else:
            if not self.paint_rect.adjusted(-DOT_RADIUS, -DOT_RADIUS,
                    DOT_RADIUS, DOT_RADIUS).contains(event.pos()):
                new_point = -1
            else:
                new_point = self.findPoint(event.pos())
            if new_point != self.hover_point:
                self.hover_point = new_point
                self.repaint()

    def wheelEvent(self, event):
        """
        Handles mouse wheel event to change data and alpha range and to
        control vertical position of color points.

        """
        try:
            delta = event.delta()
        except AttributeError:
            # Qt5
            delta = event.angleDelta().y()

        delta /= -1000.0

        for key, rect in self.text_boxes.items():
            if rect.contains(event.pos()):
                vrange = self.vmax - self.vmin
                if key == "amax":
                    self.amax = max(min(self.amax * (1.0 + delta), 1.0), 0.0)
                elif key == "vmin":
                    self.vmin = min(self.vmin + vrange * delta,
                                    self.vmax - EPS)
                else:
                    self.vmax = max(self.vmax + vrange * delta,
                                    self.vmin + EPS)
                self.repaint()
                return

        for i, p in enumerate(self.points):
            x, y, r, g, b = p
            y -= y * delta
            y = min(max(y, 0.0), 1.0)
            self.points[i] = (x, y, r, g, b)

        self.repaint()
        self.updateVolumeColors()

    def movePoints(self, event):
        """
        Moves selected point(s).
        """
        # delta == 2 if moving three points
        delta = 2 if event.modifiers() == Qt.ControlModifier else 1

        num_points = len(self.points)
        x, y, r, g, b = self.points[self.point]
        new_x = self.xToData(self.convertX(event.x()))
        dx = new_x - x

        if dx < 0:
            min_x = self.points[self.point - delta][
                0] if self.point > delta - 1 else self.vmin
            x0 = self.points[self.point - delta + 1][
                0] if self.point > 0 else x
            if x0 + dx < min_x:
                dx = min_x - x0
        else:
            max_x = self.points[self.point + delta][
                0] if self.point < num_points - delta else self.vmax
            x0 = self.points[self.point + delta - 1][
                0] if self.point < num_points - delta + 1 else x
            if x0 + dx > max_x:
                dx = max_x - x0

        new_x = x + dx
        new_y = self.yToAlpha(self.convertY(event.y()))

        # apply constrained motion
        new_x = x if self.constraint == 'y' else new_x
        new_y = y if self.constraint == 'x' or delta > 1 else new_y

        self.points[self.point] = (new_x, new_y, r, g, b)

        if delta > 1:
            dx = new_x - x
            if self.point > 0:
                x, y, r, g, b = self.points[self.point - 1]
                self.points[self.point - 1] = (x + dx, y, r, g, b)
            if self.point < num_points - 1:
                x, y, r, g, b = self.points[self.point + 1]
                self.points[self.point + 1] = (x + dx, y, r, g, b)

        self.repaint()
        s = "value: %.3f\nalpha: %.3f" % (new_x, new_y)
        QtWidgets.QToolTip.showText(self.mapToGlobal(event.pos()), s)
        if self.real_time:
            self.updateVolumeColors()

    def addPoint(self, pos, three_points):
        """
        Add a new color point.
        @param pos: Position in widget's coordinates.
        @type pos: L{QPoint}
        """
        # get color of new point
        r, g, b = next(self.color_cycle)
        new_x = self.convertX(pos.x())
        new_y = self.yToAlpha(self.convertY(pos.y()))
        new_index = len(self.points)
        for index, point in enumerate(self.points):
            if new_x < self.dataToX(point[0]):
                new_index = index
                break

        self.points.insert(new_index, (self.xToData(new_x), new_y, r, g, b))

        if three_points:
            new_x = self.convertX(pos.x()-10)
            self.points.insert(new_index, (self.xToData(new_x), 0, r, g, b))
            new_x = self.convertX(pos.x()+10)
            self.points.insert(new_index+2, (self.xToData(new_x), 0, r, g, b))
            new_index += 1

        self.point = new_index
        self.repaint()
        self.updateVolumeColors()

    def removePoints(self, three_points):
        """
        Removes one or three color points.
        """
        if self.point >= 0:
            del self.points[self.point]
            if three_points:
                if self.point < len(self.points):
                    del self.points[self.point]
                if self.point > 0:
                    del self.points[self.point - 1]
            self.point = -1
            self.repaint()
            self.updateVolumeColors()

    def findPoint(self, pos):
        """
        Finds a color point at a given cursor position.
        """
        for index, point in enumerate(self.points):
            x, y, r, g, b = point
            dx = self.dataToX(x) * (self.width() - self.left_margin
                                    ) - pos.x() + self.left_margin
            dy = (self.height() - self.bottom_margin) * (1.0 - self.alphaToY(y)
                                                         ) - pos.y()
            if dx * dx + dy * dy < 4 * DOT_RADIUS * DOT_RADIUS:
                return index
        return -1

    def setHistogram(self, histogram):
        """
        Sets histogram of the volume.

        @param histogram: Volume histogram: min value, max value, avg value,
        std deviation, value0, value1, value2...
        @type histogram: list of floats
        """
        self.vmin = histogram[0]
        self.vmax = histogram[1]

        # check for flat data (useless, but avoid errors in GUI)
        if self.vmin == self.vmax:
            self.vmin -= 1.0
            self.vmax += 1.0
        elif math.isnan(self.vmin) or math.isnan(self.vmax):
            print('Warning: setHistogram vmin={} vmax={}'.format(
                self.vmin, self.vmax))
            self.vmin = self.original_vmin
            self.vmax = self.original_vmax
            return

        self.original_vmin = self.vmin
        self.original_vmax = self.vmax
        self.path = []

        hist = histogram[4:]
        N = len(hist)
        if N == 0:
            return

        # cut extreme peaks in distribution
        shist = sorted(hist)
        q90 = shist[int(N * 0.9)]
        max_value = min(q90 * 4, shist[N - 1])
        if max_value == 0.0:
            return

        xstep = 1.0 / N
        ynorm = 1.0 / max_value
        x = 0.0
        for v in hist:
            x += xstep
            y = v * ynorm
            self.path.append((x, y))

    def setColors(self, colors):
        """
        Sets color map in the widget.

        @param colors: Flat list of values and corresponding RGBA colors
        listed as value0, R0, G0, B0, A0, value1, R1, G1, B1, A1, ...
        @type colors: list of floats
        """
        if self.ignore_set_colors:
            return
        self.points = []
        for p in range(0, len(colors), 5):
            v = colors[p]
            r = colors[p + 1]
            g = colors[p + 2]
            b = colors[p + 3]
            a = colors[p + 4]
            x = v
            y = a

            if math.isnan(x):
                print('Warning: setColors x={}'.format(x))
                return

            self.points.append((x, y, r, g, b))

        self.update()

    def reset(self):
        """
        Resets vmin and vmax to original values.
        """
        self.vmin = self.original_vmin
        self.vmax = self.original_vmax
        self.repaint()
        self.updateVolumeColors()

    def displayTextDialog(self, text):
        """
        Opens a generic text dialog and displays provided text.
        """
        if not self.text_dialog:
            self.text_dialog = QtWidgets.QDialog()
            layout = QtWidgets.QVBoxLayout()
            self.text_dialog.setLayout(layout)
            self.text_dialog.text_display = QtWidgets.QPlainTextEdit()
            self.text_dialog.text_display.setReadOnly(True)
            self.text_dialog.text_display.setStyleSheet(
                "font-family: monospace, courier")
            size = self.text_dialog.size()
            size.setWidth(DEFAULT_TEXT_DIALOG_WIDTH)
            self.text_dialog.resize(size)
            layout.addWidget(self.text_dialog.text_display)

        self.text_dialog.text_display.setPlainText(text)
        self.text_dialog.show()
        self.text_dialog.raise_()

    def displayHelp(self):
        """
        Displays text dialog with help.
        """
        self.displayTextDialog(VOLUME_HELP)

    def displayScript(self):
        """
        Displays contents of volume color ramp as script.
        """
        import random

        r = self.getColors()
        rname = 'ramp%03d' % random.randint(0, 999)
        s = ['### cut below here and paste into script ###\n']
        s.append('cmd.volume_ramp_new(%s, [\\\n' % repr(rname))
        for i in range(0, len(r), 5):
            s.append('    %6.2f, %.2f, %.2f, %.2f, %.2f, \\\n' %
                     tuple(r[i:i + 5]))
        s.append('    ])\n')
        s.append('### cut above here and paste into script ###\n')

        s += [
            '\n',
            'Paste into a .pml or .py script or your pymolrc file and use this\n',
            'named color ramp on the PyMOL command line like this:\n',
            '\n',
            'PyMOL> volume_color yourvolume, %s\n' % rname,
        ]

        self.displayTextDialog(''.join(s))

    def windowTopLevelChanged(self, floating):
        """
        Update widget colors based on floating state.
        """
        if floating:
            self.update_cb.setStyleSheet("color: black;")
        else:
            self.update_cb.setStyleSheet("color: lightgrey;")


def VolumePanelDocked(parent, *args, **kwargs):
    widget = QtWidgets.QWidget(parent)
    window = QtWidgets.QDockWidget(parent)
    _VolumePanel(widget, window, *args, **kwargs)
    window.setWidget(widget)
    parent.addDockWidget(Qt.BottomDockWidgetArea, window)

    # disabled: always use default style
    # window.topLevelChanged.connect(widget.editor.windowTopLevelChanged)

    return window


def VolumePanelDialog(parent, *args, **kwargs):
    window = QtWidgets.QDialog(parent)
    _VolumePanel(window, window, *args, **kwargs)
    return window


VolumePanel = VolumePanelDocked


class _VolumePanel(object):
    def __init__(self, widget, window=None, name='', _self=None):

        if window:
            window.setWindowTitle(name + ' - Volume Color Map Editor')

        cmd = _self

        layout = QtWidgets.QVBoxLayout()
        widget.setLayout(layout)

        widget.editor = VolumeEditorWidget(widget, name, _self)
        layout.addWidget(widget.editor)
        layout.setContentsMargins(5, 5, 5, 5)
        get_colors_btn = QtWidgets.QPushButton("Get colors as script")
        get_colors_btn.setAutoDefault(False)
        get_colors_btn.clicked.connect(widget.editor.displayScript)
        help_btn = QtWidgets.QPushButton("Help")
        help_btn.setAutoDefault(False)
        help_btn.clicked.connect(widget.editor.displayHelp)
        reset_btn = QtWidgets.QPushButton("Reset Data Range")
        reset_btn.setAutoDefault(False)
        reset_btn.clicked.connect(widget.editor.reset)
        widget.editor.update_cb = QtWidgets.QCheckBox(
            "Update volume colors in real-time")
        widget.editor.update_cb.setObjectName("volume_checkbox")
        widget.editor.update_cb.setChecked(True)
        widget.editor.update_cb.setMinimumWidth(30)
        widget.editor.update_cb.toggled.connect(
            widget.editor.toggleRealTimeUpdates)

        button_layout = QtWidgets.QHBoxLayout()
        button_layout.addWidget(get_colors_btn)
        button_layout.addWidget(reset_btn)
        button_layout.addWidget(help_btn)
        button_layout.addStretch()
        button_layout.addWidget(widget.editor.update_cb)

        layout.addLayout(button_layout)

        histogram = cmd.get_volume_histogram(name)
        widget.editor.setHistogram(histogram)

        colors = cmd.volume_color(name)
        widget.editor.setColors(colors)
