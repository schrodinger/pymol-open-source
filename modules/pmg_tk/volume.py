'''
Volume color ramp GUI
'''

import os
import math
import colorsys
import itertools

import tkinter as Tkinter

try:
    from pymol import cmd
except ImportError:
    cmd = None

CIRCLE_RADIUS = 4
STATE_SHIFT = 0x0001
STATE_CTRL  = 0x0004

DEFAULT_COLORS = [
    (1., 1., 0.),
    (1., 0., 0.),
    (0., 0., 1.),
    (0., 1., 0.),
    (0., 1., 1.),
    (1., 0., 1.),
]

help = '''
VOLUME PANEL HELP

--------------------------------------------------
Canvas Mouse Actions (no Point under Cursor)

  L-Click            Add point
  CTRL+L-Click       Add 3 points (isosurface)

--------------------------------------------------
Mouse Actions with Point under Cursor

  L-Click            Edit point color
  CTRL+L-Click       Edit 3 points

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

class ColorChooser(Tkinter.Toplevel):
    """
    HSV-based color chooser modal dialog with data and alpha entry fields.
    """

    def __init__(self, panel, points, title='Color Chooser'):
        '''
        @type panel: VolumePanel
        @type points: sequence of VRGBA
        '''
        self.panel = panel

        parent = panel.frame
        Tkinter.Toplevel.__init__(self, parent)
        self.transient(parent)
        self.withdraw()
        self.resizable(0, 0)

        self.title(title)
        self.wm_protocol('WM_DELETE_WINDOW', self.accept)

        self.points = points
        self._undo_data = [p.__getstate__() for p in points]

        imgName = os.path.expandvars("$PYMOL_DATA/pymol/hsv.ppm")
        img = Tkinter.PhotoImage(file=imgName)

        self.image_width = img.width()
        self.image_height = img.height()

        frame = Tkinter.Frame(self)
        self.canvas = Tkinter.Canvas(frame,
                width=self.image_width,
                height=self.image_height)
        self.canvas.pack(side=Tkinter.TOP)

        self.canvas.create_image((0, 0), image=img, anchor=Tkinter.NW)
        self.oval_id = self.canvas.create_oval((-10,) * 4)

        p = points[0]
        self.var_value = Tkinter.DoubleVar(self, value=p.value)
        self.var_alpha = Tkinter.DoubleVar(self, value=p.alpha)
        self.var_color = Tkinter.StringVar(self)
        self.set_rgb(p.rgb)

        self.canvas.bind("<Button-1>", self.pick_color)
        self.canvas.bind('<Button1-Motion>', self.pick_color)
        self.canvas.bind("<Double-Button-1>", self.accept)

        grid = Tkinter.Frame(frame)

        for i, text, var in [
                (0, 'Color:', self.var_color),
                (1, 'Data value:', self.var_value),
                (2, 'Opacity [0..1]:', self.var_alpha),
                ]:
            e = Tkinter.Label(grid, text=text)
            e.grid(row=i, column=0, sticky='nw')
            e = Tkinter.Entry(grid, width=7, textvariable=var)
            e.grid(row=i, column=1, sticky='nw')
            e.bind('<Return>', self.accept)
            e.bind('<FocusOut>', self.update_canvas)

        grid.pack(side=Tkinter.LEFT)

        buttonbox = Tkinter.Frame(frame)
        Tkinter.Button(buttonbox, text='OK', command=self.accept).pack(side=Tkinter.RIGHT)
        Tkinter.Button(buttonbox, text='Cancel', command=self.cancel).pack(side=Tkinter.RIGHT)

        buttonbox.pack(side=Tkinter.BOTTOM, fill='x')

        x, y = self.winfo_pointerxy()
        self.geometry("+%d+%d" % (x + 20, y - 20))
        frame.pack(padx=5, pady=5)

        self.deiconify()
        self.grab_set()
        self.focus_set()
        self._img = img # keep a reference!

    def update_upstream(self):
        '''
        Update the volume panel
        '''
        p = self.points[0]
        p.value = self.var_value.get()
        p.alpha = self.var_alpha.get()
        rgb = self.get_rgb()
        for p in self.points:
            p.rgb = rgb
        self.panel.ramp_changed()

    def move_oval(self, X, Y):
        self.canvas.coords(self.oval_id, (
            X - CIRCLE_RADIUS, Y - CIRCLE_RADIUS,
            X + CIRCLE_RADIUS, Y + CIRCLE_RADIUS))

    def update_canvas(self, event):
        '''
        Update canvas from entry fields
        '''
        self.set_rgb(self.get_rgb())

    def set_rgb(self, rgb, update_oval=True):
        '''
        Update RGB entry field and optionally move the oval
        '''
        self.var_color.set('#%02x%02x%02x' % tuple(
            int(min(max(v, 0), 1) * 0xFF) for v in rgb))

        if update_oval:
            h, s, v = colorsys.rgb_to_hsv(*rgb)
            X = h * self.image_width
            H2 = self.image_height / 2
            Y = s * H2 if (s < 1.0) else (2.0 - v) * H2
            self.move_oval(X, Y)

        if self.panel.instant_update.get():
            self.update_upstream()

    def get_rgb(self):
        '''
        return RGB float (0..1) tuple
        '''
        s = self.var_color.get()
        return tuple(int(s[i:i+2], 16) / 255. for i in range(1, 6, 2))

    def pick_color(self, event):
        '''
        Move the oval and update the RGB entry field
        '''
        x = min(max(event.x, 0), self.image_width)
        y = min(max(event.y, 0), self.image_height)
        self.move_oval(x, y)
        h = float(x) / self.image_width
        s, v = 2.0 * y / self.image_height, 1.0
        if s > 1.0:
            s, v = 1.0, 2.0 - s
        self.set_rgb(colorsys.hsv_to_rgb(h, s, v), False)

    def cancel(self, event=None):
        for p, undo in zip(self.points, self._undo_data):
            p.__setstate__(undo)
        self.panel.ramp_changed()
        self.destroy()

    def accept(self, event=None):
        self.update_upstream()
        self.destroy()

class VRGBA(object):
    '''
    Simple Value-RGB-Alpha type
    '''
    def __init__(self, canvas, value, alpha, rgb):
        '''
        @type canvas: Tkinter.Canvas
        @type value: float
        @type alpha: float
        @type rgb: tuple
        @param rgb: 0..1 float RGB
        '''
        self.value = value
        self.rgb = rgb
        self.alpha = alpha
        self.canvas_id = canvas.create_oval((-10,) * 4, tags=('dynamic', 'colorpoint'))

    def __iter__(self):
        yield self.value
        yield self.rgb[0]
        yield self.rgb[1]
        yield self.rgb[2]
        yield self.alpha

    def __getstate__(self):
        return list(self)

    def __setstate__(self, state):
        self.value = state[0]
        self.rgb = state[1:4]
        self.alpha = state[4]

    @property
    def alpha(self):
        '''
        @type: float
        '''
        return self._alpha

    @alpha.setter
    def alpha(self, v):
        self._alpha = min(max(v, 0.0), 1.0)

    @property
    def hexcolor(self):
        '''
        @type: str
        '''
        return "#%02x%02x%02x" % tuple(int(v * 0xFF) for v in self.rgb)

class RangeEntry(Tkinter.Entry):
    '''
    Entry field for canvas value range
    '''
    def __init__(self, parent, panel, vname, rname):
        '''
        @type parent: Tkinter.Widget
        @type panel: VolumePanel
        @type vname: str
        @param vname: panel field name of linked property
        @type rname: str
        @param vname: panel field name of corresponding range
        '''
        Tkinter.Entry.__init__(self, parent, width=5, bg='white')
        self._panel = panel
        self._vname = vname
        self._rname = rname
        self.bind('<Return>', self.onchange)
        self.bind('<FocusOut>', self.onchange)
        self.bind('<Button-4>', self.increment)
        self.bind('<Button-5>', self.decrement)
        self.bind('<Up>', self.increment)
        self.bind('<Down>', self.decrement)
        self._update()

    def _set(self, v):
        setattr(self._panel, self._vname, v)
        self._update()
        self._panel.redraw()

    def _update(self):
        v = getattr(self._panel, self._vname)
        self.delete(0, Tkinter.END)
        self.insert(0, v)
        self.config(bg='white')

    def increment(self, event, m=0.05):
        r = getattr(self._panel, self._rname)
        v = getattr(self._panel, self._vname)
        self._set(v + r * m)
        return 'break'

    def decrement(self, event, m=0.05):
        return self.increment(event, -m)

    def onchange(self, event):
        try:
            v_new = float(self.get())
        except ValueError:
            self.config(bg='red')
            return False
        v_old = getattr(self._panel, self._vname)
        if v_new != v_old:
            self._set(v_new)
        return True

class VRGBACanvas(object):
    '''
    Editable Value-RGBA ramp graph
    '''

    def __init__(self, parent):
        '''
        @type parent: Tkinter.Widget
        @type name: name of volume in PyMOL
        '''
        self.frame = Tkinter.Frame(parent)
        self.canvas = canvas = Tkinter.Canvas(self.frame)
        self.buttonbox = Tkinter.Frame(self.frame)

        self.canvas_ids = {} # ID -> VRGBA

        self.ramp = []
        self.instant_update = Tkinter.BooleanVar(parent, value=1)

        self.canvas_width = 0
        self.canvas_height = 0

        self._vmin = -5.
        self._vmax = 5.
        self._amin = 0.
        self._amax = 1.

        self._histcoords = ()
        self.drag_ids = None

        self.color_cycle = itertools.cycle(DEFAULT_COLORS)

        # plot area padding
        self.pad_left = 50
        self.pad_right = 10
        self.pad_bottom = 30
        self.pad_top = 10

        # entries
        self.minmax_entries = [
            RangeEntry(self.canvas, self, 'vmin', 'vrange'),
            RangeEntry(self.canvas, self, 'vmax', 'vrange'),
            RangeEntry(self.canvas, self, 'amax', 'arange'),
        ]
        self.minmax_ids = [
            canvas.create_window((-10, -10), window=self.minmax_entries[i], anchor=a)
            for i, a in enumerate(['nw', 'ne', 'ne'])
        ]

        # color line
        self.colorline_id = canvas.create_line((-10,) * 4, tags=('dynamic',))

        # coord tooltip
        self.tooltip_id = canvas.create_text([-1, -1], text='', tags=('dynamic',), anchor='sw')

        # bindings
        canvas.bind('<Configure>', self.onconfig)
        for i in range(1, 4):
            self.canvas.bind('<Button-%d>' % i, self.onmousedown)
            self.canvas.bind('<ButtonRelease-%d>' % i, self.onmouseup)
            self.canvas.bind('<Button%d-Motion>' % i, self.ondrag)
        for i in range(4, 6):
            self.canvas.bind('<Button-%d>' % i, self.onmousewheel)

        # buttons
        Tkinter.Button(self.buttonbox, text='Get colors as script',
                command=self.popup_script).pack(side="left")
        Tkinter.Button(self.buttonbox, text='Help',
                command=self.popup_help).pack(side="left")
        Tkinter.Checkbutton(self.buttonbox, text='Update Volume while dragging',
                var=self.instant_update).pack(side="right")

        self.buttonbox.pack(side='bottom', fill='x')
        self.canvas.pack(side='left', expand=1, fill='both')

    def popup_help(self, event=None):
        text_dialog(self.frame, help, 'Volume Panel Help')

    def popup_script(self, event=None):
        '''
        Open a dialog with a script equivalent of this color ramp
        '''
        import random

        r = self.get_flat()
        rname = 'ramp%03d' % random.randint(0, 999)
        s = ['### cut below here and paste into script ###']
        s.append('cmd.volume_ramp_new(%s, [\\' % repr(rname))
        for i in range(0, len(r), 5):
            s.append('    %6.2f, %.2f, %.2f, %.2f, %.2f, \\' % tuple(r[i:i+5]))
        s.append('    ])')
        s.append('### cut above here and paste into script ###')

        s += [
            '',
            'Paste into a .pml or .py script or your pymolrc file and use this',
            'named color ramp on the PyMOL command line like this:',
            '',
            'PyMOL> volume_color yourvolume, %s' % rname,
        ]

        text_dialog(self.frame, '\n'.join(s), 'Script')

    def onhover(self, event=None):
        '''
        Mouse enter/leave event handler. Changes the mouse cursor when hovering over color points.
        '''
        if event and event.type == '7':
            cursor = 'crosshair'
        else:
            cursor = ''
        self.canvas.config(cursor=cursor)

    def bind_hover_events(self):
        '''
        Bind mouse enter/leave events for color points
        '''
        self.canvas.tag_bind('colorpoint', '<Enter>', self.onhover, add=False)
        self.canvas.tag_bind('colorpoint', '<Leave>', self.onhover, add=False)

    def onconfig(self, event):
        '''
        Event callback for canvas configure, redraws everything
        '''
        self.canvas_width = event.width - self.pad_left - self.pad_right
        self.canvas_height = event.height - self.pad_top - self.pad_bottom
        self.redraw()

    def set_flat(self, ramp):
        '''
        Set color points from flat (v0, r0, g0, b0, a0, v1, r1, ...) list
        @type ramp: sequence of floats
        '''
        self.canvas.delete('colorpoint')
        self.ramp = [VRGBA(self.canvas, ramp[i], ramp[i + 4], ramp[i + 1:i + 4])
                for i in range(0, len(ramp), 5)]
        self.update_canvas_ids()
        self.plot_ramp()

    def update_canvas_ids(self):
        '''
        Update mapping canvas_ids -> color point
        '''
        self.canvas_ids = dict((p.canvas_id, p) for p in self.ramp)
        self.bind_hover_events()

    def get_flat(self):
        '''
        Get color points as flat (v0, r0, g0, b0, a0, v1, r1, ...) list
        @rtype: list of float
        '''
        flat = []
        for p in self.ramp:
            flat.append(p.value)
            flat.extend(p.rgb)
            flat.append(p.alpha)
        return flat

    @property
    def vmin(self):
        return self._vmin

    @property
    def vmax(self):
        return self._vmax

    @vmin.setter
    def vmin(self, v):
        self._vmin = min(v, self._vmax - 1e-2)

    @vmax.setter
    def vmax(self, v):
        self._vmax = max(v, self._vmin + 1e-2)

    @property
    def vrange(self):
        '''
        value range of ploting area
        @type: float
        '''
        return self.vmax - self.vmin

    @property
    def amin(self):
        '''
        alpha minimum of plotting area
        @type: float
        '''
        return self._amin

    @property
    def amax(self):
        '''
        alpha maximum of plotting area
        @type: float
        '''
        return self._amax

    @amax.setter
    def amax(self, a):
        self._amax = min(max(a, 0.1), 1.0)

    @property
    def arange(self):
        '''
        alpha range of plotting area
        @type: float
        '''
        return self.amax - self.amin

    def add_point(self, value, alpha, color):
        '''
        Add a RGBA color at given value
        @type value: float
        @type alpha: float
        @type color: 3-tuple
        '''
        pt = VRGBA(self.canvas, value, alpha, color)

        for i, c in enumerate(self.ramp):
            if c.value > value:
                self.ramp.insert(i, pt)
                break
        else:
            self.ramp.append(pt)
        self.update_canvas_ids()

    def closest_point(self, value):
        '''
        Search for point closest to value
        @type value: float
        @rtype: VRGBA or None
        '''
        if not self.ramp:
            return None
        return min(self.ramp, key=lambda p: abs(p.value - value))

    def remove_point(self, p):
        '''
        @type p: VRGBA
        '''
        self.canvas.delete(p.canvas_id)
        self.ramp.remove(p)
        self.update_canvas_ids()
        self.onhover()

    # coordinate transform methods

    def valueToCanvas(self, v):
        x = (v - self.vmin) / self.vrange
        return self.pad_left + self.canvas_width * x

    def canvasToValue(self, x):
        x = float(x - self.pad_left) / self.canvas_width
        return self.vrange * x + self.vmin

    def alphaToCanvas(self, a):
        a = min(max(a, 0.), 1.)
        y = (a - self.amin) / self.arange
        y = math.log(1.0 + 9.0 * y, 10.0)
        return self.pad_top + self.canvas_height * (1.0 - y)

    def canvasToAlpha(self, y):
        y = 1.0 - float(y - self.pad_top) / self.canvas_height
        y = (10.0 ** y - 1.0) / 9.0
        a = y * self.arange + self.amin
        return min(max(a, 0.), 1.)

    def vaToCanvas(self, coords):
        return self._transCoords(coords,
                self.valueToCanvas, self.alphaToCanvas)

    def canvasToVA(self, coords):
        return self._transCoords(coords,
                self.canvasToValue, self.canvasToAlpha)

    def _transCoords(self, coords, xfunc, yfunc):
        c = [0] * len(coords)
        c[0::2] = [xfunc(x) for x in coords[0::2]]
        c[1::2] = [yfunc(y) for y in coords[1::2]]
        return c

    # end coordinate transformation methods

    def create_line_va(self, coords, *args, **kwargs):
        '''
        Draw a line with value and alpha coordinates
        @type coords: sequence of floats
        '''
        c = self.vaToCanvas(coords)
        self.canvas.create_line(c, *args, **kwargs)

    def set_hist(self, hist):
        '''
        Set the histogram, see also C{plot_hist}.
        @type coords: sequence of floats
        @param coords: histogram with (min, max, mean, stdev) as first 4
        elements, followed by equally spaced bin counts.
        '''
        vmin, vmax = hist[:2]
        hist = hist[4:]
        try:
            ihistmax = 1.0 / max(hist)
        except ZeroDivisionError:
            ihistmax = 0.0
        binwidth = (vmax - vmin) / float(len(hist) - 1)
        if not self._histcoords:
            self.vmin = vmin
            self.vmax = vmax
            for e in self.minmax_entries:
                e._update()
        c = self._histcoords = []
        for i in range(len(hist)):
            c.append(float(i) * binwidth + vmin)
            c.append(hist[i] * ihistmax)

    def plot_hist(self):
        '''
        Plot the histogram which was set by C{set_hist}
        '''
        self.canvas.delete('hist')
        self.create_line_va(self._histcoords, fill='blue', tags=('hist'))

    def plot_axis(self):
        '''
        Plot the coordinate system axes
        '''
        self.canvas.delete('axis')

        width, height = self.canvas_width, self.canvas_height
        if width < 1 or height < 1:
            return

        x0 = self.valueToCanvas(self.vmin)
        x1 = self.valueToCanvas(self.vmax)
        y0 = self.alphaToCanvas(self.amin)
        y1 = self.alphaToCanvas(self.amax)
        gridkw = {'dash': (2, 5), 'fill': '#999999', 'tags': ('axis')}

        # x-grid
        u30 = self.canvasToValue(30 + self.pad_left) - self.vmin
        try:
            u = 5 ** int(1 + math.log(u30, 5))
        except ValueError:
            u = self.vrange / 4.0
        v = round(self.vmin / u) * u
        while v < self.vmax:
            v += u
            x = self.valueToCanvas(v)
            self.canvas.create_line([x, y0, x, y1], **gridkw)
            if x > self.pad_left + 50 and x < width:
                self.canvas.create_line([x, y0, x, y0 + 5], tags=('axis'))
                self.canvas.create_text([x, y0 + 7], text='%.4G' % v, tags=('axis'), anchor='n')

        # y-grid
        for i in range(1, 11):
            a = i * 0.1
            y = self.alphaToCanvas(a)
            self.canvas.create_line([x0, y, x1, y], **gridkw)
            if y > 50:
                self.canvas.create_line([x0, y, x0 - 5, y], tags=('axis'))
                self.canvas.create_text([x0 - 7, y], text='%.1f' % a, tags=('axis'), anchor='e')

        # axes
        coords = [x0, y1, x0, y0, x1, y0]
        self.canvas.create_line(coords, width=3, tags=('axis'))

        # entries
        for xy, e in [
                ((x0, y0 + 4), self.minmax_ids[0]),
                ((x1, y0 + 4), self.minmax_ids[1]),
                ((x0 - 2, y1), self.minmax_ids[2]),
                ]:
            self.canvas.coords(e, xy)

    def redraw(self):
        '''
        Redraw or reposition all elements on canvas
        '''
        self.plot_hist()
        self.plot_axis()
        self.plot_ramp()
        self.canvas.tag_raise('dynamic')

    def plot_ramp(self):
        '''
        Update color points and line coordinates
        '''
        radius = CIRCLE_RADIUS

        if self.canvas_width <= 0 or self.canvas_height <= 0:
            return

        polycoords = []

        for idx, p in enumerate(self.ramp):
            x = self.valueToCanvas(p.value)
            y = self.alphaToCanvas(p.alpha)

            polycoords.extend((x, y))

            # plot the pt
            self.canvas_set_position(p.canvas_id, x, y)
            self.canvas.coords(p.canvas_id, (x - radius, y - radius, x + radius, y + radius))
            self.canvas.itemconfig(p.canvas_id, fill=p.hexcolor)

        if len(polycoords) < 4:
            polycoords = (-10,) * 4
        self.canvas.coords(self.colorline_id, tuple(polycoords))

    def plot_colorline(self):
        '''
        Update color line coordinates
        '''
        if len(self.ramp) < 2:
            polycoords = (-10,) * 4
        else:
            polycoords = []
            for p in self.ramp:
                xy = self.canvas_get_position(p.canvas_id)
                polycoords.extend(xy)
        self.canvas.coords(self.colorline_id, tuple(polycoords))

    def canvas_get_position(self, ID):
        '''
        Get x, y canvas position (center) of color point
        @type ID: int
        '''
        radius = CIRCLE_RADIUS
        coords = self.canvas.coords(ID)
        return coords[0] + radius, coords[1] + radius

    def canvas_set_position(self, ID, x, y):
        '''
        Set canvas position (center) of color point
        '''
        x_old, y_old = self.canvas_get_position(ID)
        self.canvas.move(ID, x - x_old, y - y_old)

    def update_va_from_canvas(self, IDs):
        '''
        Update color points value and alpha from positions on canvas
        @type IDs: iterable
        @param IDs: canvas element ids
        '''
        for ID in IDs:
            p = self.canvas_ids[ID]
            x, y = self.canvas_get_position(ID)
            p.value = self.canvasToValue(x)
            p.alpha = self.canvasToAlpha(y)

    def canvas_remove_point(self, ID):
        '''
        @type ID: int
        @param ID: canvas element id
        '''
        p = self.canvas_ids[ID]
        return self.remove_point(p)

    def canvas_add_point(self, x, y):
        '''
        Add point on canvas. Use the color of the closest neighbor.
        @type x: int
        @type y: int
        '''
        v = self.canvasToValue(x)
        a = self.canvasToAlpha(y)
        c = next(self.color_cycle)
        self.add_point(v, a, c)

    def canvas_add_peak(self, x, y):
        '''
        Add three points on canvas, one at (x,y) and two at (x +/- 10, y=0)
        @type x: int
        @type y: int
        '''
        v1 = self.canvasToValue(x - 10)
        v2 = self.canvasToValue(x)
        v3 = self.canvasToValue(x + 10)
        a = self.canvasToAlpha(y)
        c = next(self.color_cycle)
        self.add_point(v1, 0, c)
        self.add_point(v2, a, c)
        self.add_point(v3, 0, c)

    def canvas_pick(self, x, y, pad=2.):
        '''
        Get canvas ID of color point at (x, y)
        @type x: int
        @type y: int
        @type pad: float
        @param pad: extra padding for easier picking
        '''
        a = self.canvas.find_closest(x, y)
        if not a or a[0] not in self.canvas_ids:
            return None
        center = self.canvas_get_position(a[0])
        d2 = (center[0] - x) ** 2 + (center[1] - y) ** 2
        if d2 > (CIRCLE_RADIUS + pad) ** 2:
            return None
        return a[0]

    def onmousewheel(self, event):
        '''
        Mouse wheel event handler
        '''
        factor = 1.05 if event.num == 4 else 0.95;

        for p in self.ramp:
            p.alpha *= factor

        self.ramp_changed()

    def onmousedown(self, event):
        '''
        Mouse down event handler
        '''
        x, y = event.x, event.y

        canvas_right = self.pad_left + self.canvas_width

        self.dragging = False
        self.drag_ids = ()
        self.drag_button = event.num
        self.dragging_x = self.drag_start_x = x
        self.dragging_y = self.drag_start_y = y
        c = self.canvas_pick(x, y)

        if not c:
            if not event.num == 1:
                return
            # add
            if x < self.pad_left \
                    or x > canvas_right \
                    or y < self.pad_top \
                    or y > self.pad_top + self.canvas_height:
                return

            if event.state & STATE_CTRL:
                self.canvas_add_peak(event.x, event.y)
            else:
                self.canvas_add_point(event.x, event.y)
            self.ramp_changed()
        else:
            # move
            self.drag_ids = [c]
            i = self.ramp.index(self.canvas_ids[c])
            i_lower = i - 1
            i_upper = i + 1
            x, y = self.canvas.coords(c)[:2]
            dx_upper, dy = event.x - x, event.y - y - CIRCLE_RADIUS
            dx_lower = dx_upper

            if event.state & STATE_CTRL:
                # move peak
                if i > 0:
                    c = self.ramp[i - 1].canvas_id
                    self.drag_ids.append(c)
                    x = self.canvas.coords(c)[0]
                    dx_lower = event.x - x
                    i_lower = i - 2
                if i < len(self.ramp) - 1:
                    c = self.ramp[i + 1].canvas_id
                    self.drag_ids.append(c)
                    x = self.canvas.coords(c)[0]
                    dx_upper = event.x - x
                    i_upper = i + 2
                self.allowed_y = (event.y, event.y)
            else:
                # move single point
                self.allowed_y = (
                        min(event.y, self.pad_top + dy),
                        max(event.y, self.pad_top + dy + self.canvas_height))

            # set up allowed x range
            x_lower = x_upper = None
            if i_lower > -1:
                x = self.canvas.coords(self.ramp[i_lower].canvas_id)[0]
                x_lower = x + dx_lower
            else:
                x_lower = self.pad_left + dx_lower - CIRCLE_RADIUS
            if i_upper < len(self.ramp):
                x = self.canvas.coords(self.ramp[i_upper].canvas_id)[0]
                x_upper = x + dx_upper
            else:
                x_upper = canvas_right + dx_lower - CIRCLE_RADIUS
            self.allowed_x = (min(event.x, x_lower), max(event.x, x_upper))

        # middle-click or SHIFT-click remove
        if event.num == 2 or event.num == 1 and event.state & STATE_SHIFT:
            for c in self.drag_ids:
                self.canvas_remove_point(c)
                self.ramp_changed()
            self.drag_ids = ()

    def onmouseup(self, event):
        '''
        Mouse up event handler
        '''
        if not self.drag_ids:
            return

        if self.dragging:
            self.canvas.itemconfig(self.tooltip_id, text='')
            if not self.instant_update.get():
                self.update_va_from_canvas(self.drag_ids)
        else:
            points = [self.canvas_ids[ID] for ID in self.drag_ids]
            chooser = ColorChooser(self, points)

        self.drag_ids = ()

    def ramp_changed(self):
        self.plot_ramp()

    def clamp_xy_allowed(self, x, y):
        '''
        Clamp x and y to allowed dragging range

        If dragging with right mouse button, drag in x or y direction
        only (choose axis closer to pointer).
        '''
        if self.allowed_x[0] is not None and self.allowed_x[0] > x:
            x = self.allowed_x[0]
        if self.allowed_x[1] is not None and self.allowed_x[1] < x:
            x = self.allowed_x[1]
        if self.allowed_y[0] is not None and self.allowed_y[0] > y:
            y = self.allowed_y[0]
        if self.allowed_y[1] is not None and self.allowed_y[1] < y:
            y = self.allowed_y[1]
        if self.drag_button == 3:
            if abs(x - self.drag_start_x) < abs(y - self.drag_start_y):
                x = self.drag_start_x
            else:
                y = self.drag_start_y
        return x, y

    def ondrag(self, event):
        '''
        Mouse dragging event handler
        '''
        if not self.drag_ids:
            return
        self.dragging = True

        x, y = self.clamp_xy_allowed(event.x, event.y)

        dx = x - self.dragging_x
        dy = y - self.dragging_y
        self.dragging_x = x
        self.dragging_y = y

        # move points
        for c in self.drag_ids:
            self.canvas.move(c, dx, dy)

        # tooltip
        x, y = self.canvas_get_position(self.drag_ids[0])
        v = self.canvasToValue(x)
        a = self.canvasToAlpha(y)
        dx = dy = CIRCLE_RADIUS
        if y < 40:
            anchor = 'n'
        else:
            anchor = 's'
            dy *= -1
        if x > self.pad_left + self.canvas_width - 80:
            anchor += 'e'
            dx *= -1
        else:
            anchor += 'w'
        self.canvas.itemconfig(self.tooltip_id, anchor=anchor,
                text='value: %G\nalpha: %.2f' % (v, a))
        self.canvas.coords(self.tooltip_id, (x + dx, y + dy))

        # colorline
        self.plot_colorline()
 
        if self.instant_update.get():
            self.update_va_from_canvas(self.drag_ids)

    def pack(self, cnf={}, **kw):
        '''
        Tkinter: pack in parent widget
        '''
        c = {'expand': 1, 'fill': 'both'}
        c.update(cnf)
        c.update(kw)
        self.frame.pack(c)

def text_dialog(parent, text, title=''):
    '''
    Simple Text dialog
    '''
    import Pmw
    dialog = Pmw.TextDialog(parent, title=title)
    dialog.insert('end', text)

class VolumePanel(VRGBACanvas):
    '''
    Editable volume ramp graph
    '''

    def __init__(self, parent, name, _self=cmd):
        '''
        @type parent: Tkinter.Widget
        @type name: name of volume in PyMOL
        '''
        self._super = super(VolumePanel, self)
        self._super.__init__(parent)

        self.cmd = _self
        self.name = name

        # set data from volume
        hist = _self.get_volume_histogram(name)
        ramp = _self.volume_color(name)
        self.set_hist(hist)
        self.set_flat(ramp)

    def update_volume(self):
        '''
        Update volume colors in PyMOL
        '''
        ramp = self.get_flat()
        self.cmd.volume_color(self.name, ramp, _guiupdate=False)

    def update_va_from_canvas(self, *args):
        self._super.update_va_from_canvas(*args)
        self.update_volume()

    def ramp_changed(self):
        self._super.ramp_changed()
        self.update_volume()

if __name__ == '__main__':
    ramp = [-2.4, 0.0, 0.0, 1.0, 0.0,
            1.7, 0.0, 0.0, 1.0, 0.0,
            1.9, 1.0, 0.0, 0.2,
            0.2, 2.12, 0.0, 0.0, 1.0, 0.0,
            4.99, 0.0, 0.0, 1.0, 0.0]
    hist = [-20.0, 20.0, -2.63597939920146e-07, 0.9998476505279541,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0,
            80.0, 189.0, 895.0, 2571.0, 7627.0, 14798.0, 26494.0,
            32283.0, 46630.0, 51979.0, 49263.0, 42430.0, 32592.0,
            23094.0, 15904.0, 11534.0, 8759.0, 5858.0, 6157.0,
            5183.0, 4335.0, 4264.0, 3498.0, 2845.0, 2863.0, 2334.0,
            2085.0, 1719.0, 1314.0, 1047.0, 819.0, 687.0, 599.0,
            335.0, 319.0, 197.0, 142.0, 99.0, 60.0, 32.0, 32.0,
            26.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    window = Tkinter.Tk()
    canvasramp = VRGBACanvas(window)
    canvasramp.set_hist(hist)
    canvasramp.set_flat(ramp)
    canvasramp.pack()
