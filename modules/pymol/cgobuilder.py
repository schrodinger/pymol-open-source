# A* -------------------------------------------------------------------
# B* This file contains source code for the PyMOL computer program
# C* Copyright (c) Schrodinger, LLC.
# D* -------------------------------------------------------------------
# E* It is unlawful to modify or remove this copyright notice.
# F* -------------------------------------------------------------------
# G* Please see the accompanying LICENSE file for further information.
# H* -------------------------------------------------------------------
# I* Additional authors of this source file include:
# -*
# -*
# -*
# Z* -------------------------------------------------------------------

import dataclasses
import math
from typing import List, Dict

from chempy import cpv

from pymol import cgo

from . import cmd


@dataclasses.dataclass
class Point:
    """
    Represents a 3D point in space

    Attributes:
        x (float): x-coordinate
        y (float): y-coordinate
        z (float): z-coordinate
    """
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def array(self) -> List[float]:
        """
        Converts point into 3-component array
        :returns : array of xyz coordinates
        """
        return [self.x, self.y, self.z]


@dataclasses.dataclass
class Color:
    """
    Represents a Color in RGBA space

    Attributes:
        r (float): red channel value
        g (float): green channel value
        b (float): blue channel value
        a (float): alpha channel value
    """
    r: float = 0.0
    g: float = 0.0
    b: float = 0.0
    a: float = 1.0

    def array(self) -> List[float]:
        """
        Converts color into 3-component array
        """
        return [self.r, self.g, self.b]


@dataclasses.dataclass
class CGO:
    """
    Base class for all CGOs

    For all CGO types, whenever attributes are changed, `rebuild` should be
    called to update the data cache.
    """
    _data: List[float] = dataclasses.field(default_factory=list)

    def __post_init__(self):
        """
        Rebuilds CGO. All child classes should implement this.
        """
        self.rebuild()

    @property
    def data(self):
        return self._data


@dataclasses.dataclass
class Sphere(CGO):
    """
    Sphere CGO

    Attributes:
        center (Point): center position
        radius (float): sphere radius
        color (Color): sphere color
    """
    center: Point = Point(0, 0, 0)
    radius: float = 0.0
    color: Color = Point(0, 0, 0)

    def rebuild(self) -> None:
        """
        Rebuilds sphere
        """
        self._data.clear()
        self._data.append(cgo.COLOR)
        self._data.extend(self.color.array())
        self._data.append(cgo.SPHERE)
        self._data.extend(self.center.array())
        self._data.append(self.radius)


@dataclasses.dataclass
class Cylinder(CGO):
    """
    Cylinder CGO

    Attributes:
        point1 (Point): center position of first side
        point2 (Point): center position of second side
        radius (float): sphere radius
        color1 (Color): color of first cylinder half
        color2 (Color): color of second cylinder half
    """
    point1: Point = Point(0, 0, 0)
    point2: Point = Point(1, 1, 1)
    radius: float = 1.0
    color1: Color = Point(0, 0, 0)
    color2: Color = Color(0, 0, 0)

    def rebuild(self) -> None:
        """
        Rebuilds cylinder
        """
        self._data.clear()
        self._data.append(cgo.CYLINDER)
        self._data.extend(self.point1.array())
        self._data.extend(self.point2.array())
        self._data.append(self.radius)
        self._data.extend(self.color1.array())
        self._data.extend(self.color2.array())


@dataclasses.dataclass
class Torus(CGO):
    """
    Torus CGO

    Attributes:
        center (Point): torus center position
        normal (Point): torus normal direction
        radius (float): ring radius
        color (Color): torus color
        cradius (float): torus ring radius
        samples (int): number of samples to generate torus segment
        csamples (int): number of samples to generate torus ring
    """
    center: Point = Point(0.0, 0.0, 0.0)
    normal: Point = Point(0.0, 0.0, 1.0)
    radius: float = 1.0
    color: Color = Color(0.0, 1.0, 1.0, 1.0)
    cradius: float = 0.25
    samples: int = 20
    csamples: int = 20

    def rebuild(self) -> None:
        """
        Rebuilds torus
        """
        obj = []

        axis = cpv.cross_product(self.normal.array(), (0., 0., 1.))
        angle = -cpv.get_angle(self.normal.array(), (0., 0., 1.))
        matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

        def obj_vertex(x, y, z):
            return [cgo.VERTEX] + cpv.add(self.center.array(),
                                          cpv.transform(matrix, [x, y, z]))

        def obj_normal(x, y, z):
            return [cgo.NORMAL] + cpv.transform(matrix, [x, y, z])

        r = self.radius
        cr = self.cradius
        rr = 1.5 * cr
        dv = 2 * math.pi / self.csamples
        dw = 2 * math.pi / self.samples
        v = 0.0
        w = 0.0

        while w < 2 * math.pi:
            v = 0.0
            c_w = math.cos(w)
            s_w = math.sin(w)
            c_wdw = math.cos(w + dw)
            s_wdw = math.sin(w + dw)

            obj.append(cgo.BEGIN)
            obj.append(cgo.TRIANGLE_STRIP)

            obj.append(cgo.COLOR)
            obj.extend(self.color.array())

            while v < 2 * math.pi + dv:
                c_v = math.cos(v)
                s_v = math.sin(v)
                c_vdv = math.cos(v + dv)
                s_vdv = math.sin(v + dv)
                obj.extend(
                    obj_normal((r + rr * c_v) * c_w - (r + cr * c_v) * c_w,
                               (r + rr * c_v) * s_w - (r + cr * c_v) * s_w,
                               (rr * s_v - cr * s_v)))
                obj.extend(
                    obj_vertex((r + cr * c_v) * c_w, (r + cr * c_v) * s_w,
                               cr * s_v))
                obj.extend(
                    obj_normal(
                        (r + rr * c_vdv) * c_wdw - (r + cr * c_vdv) * c_wdw,
                        (r + rr * c_vdv) * s_wdw - (r + cr * c_vdv) * s_wdw,
                        rr * s_vdv - cr * s_vdv))
                obj.extend(
                    obj_vertex((r + cr * c_vdv) * c_wdw,
                               (r + cr * c_vdv) * s_wdw, cr * s_vdv))
                v += dv

            obj.append(cgo.END)
            w += dw

        self._data = obj


@dataclasses.dataclass
class Cone(CGO):
    """
    Cone CGO

    Attributes:
        tip (Point): cone tip position
        base_center (Point): center position of cone base
        radius (float): cone base radius
        color (Color): cone color
    """
    tip: Point = Point(1.0, 0.0, 0.0)
    base_center: Point = Point(0.0, 0.0, 0.0)
    radius: float = 1.0
    color: Color = Color(0.0, 0.0, 0.0)

    def rebuild(self) -> None:
        """
        Rebuilds cone
        """
        self._data = []
        self._data.append(cgo.CONE)
        self._data.extend(self.tip.array())
        self._data.extend(self.base_center.array())
        self._data.extend([self.radius, 0.0])
        self._data.extend(self.color.array())
        self._data.extend(self.color.array())
        self._data.extend([1.0, 0.0])


@dataclasses.dataclass
class Bezier(CGO):
    """
    Cubic Bezier CGO

    Attributes:
        control_pt_A (Point): on-curve bezier spline position (from)
        A_right_handle (Point): influences curve from A to B
        control_pt_B (Point): on-curve bezier spline position (to)
        B_left_handle (Point): influences curve from A to B
        A_color (Color): color at point A
        B_color (Color): color at point B
    """
    control_pt_A = Point(-5.0, 0.0, 0.0)
    A_right_handle = Point(0.0, 10.0, 0.0)
    B_left_handle = Point(1.0, -10.0, 0.0)
    control_pt_B = Point(5.0, 0.0, 0.0)

    def rebuild(self) -> None:
        """
        Rebuilds bezier spline
        """
        self._data = []
        self._data.append(cgo.BEZIER)
        self._data.extend(self.control_pt_A.array())
        self._data.extend(self.A_right_handle.array())
        self._data.extend(self.B_left_handle.array())
        self._data.extend(self.control_pt_B.array())


class CGOBuilder:
    """
    PyMOL CGO Builder

    Holds a storage of CGOs to load into PyMOL as a single object.

    Attributes:
        cgos (Dict[CGO]): Collection of CGOs
    """
    cgos: Dict[int, CGO] = {}
    _token_count = 0

    def add(self, cgo: CGO) -> int:
        """
        Adds CGO to builder

        :param cgo (CGO): cgo to add

        :return: token ID for added CGO
        """
        token = self._token_count
        self.cgos[token] = cgo
        self._token_count += 1
        return token

    def remove(self, token: int) -> None:
        """
        Removes CGO from builder

        :param token (int): token of CGO to remove
        """
        self.cgos.pop(token, None)

    def clear(self) -> None:
        """
        Clears the CGO Collection
        """
        self.cgos.clear()

    def load(self, name: str, flush: bool = True, _self=cmd) -> None:
        """
        Loads the CGO as an object into pymol

        :param name (str): Name of object
        :param flush (bool): Should reset the object data cache after loading
        """
        all_serialized = []
        for _, cgo_obj in self.cgos.items():
            all_serialized.extend(cgo_obj.data)

        _self.load_cgo(all_serialized, name)
        if flush:
            self.cgos.clear()
