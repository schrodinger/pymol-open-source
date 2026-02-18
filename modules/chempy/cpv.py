#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------
# Generic vector and matrix routines for 3-Space
# Assembled for usage in PyMOL and Chemical Python
#
# Assumes row-major matrices and arrays
# [ [vector 1], [vector 2], [vector 3] ]
#
# Raises ValueError when given bad input
#
# TODO: documentation!

from typing import Sequence
import math
import random

RSMALL4 = 0.0001

#------------------------------------------------------------------------------
def get_null() -> list[float]:
    return [0.0,0.0,0.0]

#------------------------------------------------------------------------------
def get_identity() -> list[list[float]]:
    return [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]

#------------------------------------------------------------------------------
def distance_sq(v1: Sequence[float], v2: Sequence[float]) -> float:
    d0 = v2[0] - v1[0]
    d1 = v2[1] - v1[1]
    d2 = v2[2] - v1[2]
    return (d0*d0) + (d1*d1) + (d2*d2)

#------------------------------------------------------------------------------
def distance(v1: Sequence[float], v2: Sequence[float]) -> float:
    d0 = v2[0] - v1[0]
    d1 = v2[1] - v1[1]
    d2 = v2[2] - v1[2]
    return math.sqrt((d0*d0) + (d1*d1) + (d2*d2))

#------------------------------------------------------------------------------
def length(v: Sequence[float]) -> float:
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

#------------------------------------------------------------------------------
def _random_vector() -> list[float]:
    r = random.random
    return [r()-0.5,r()-0.5,r()-0.5]

def random_displacement(v: Sequence[float], radius: float) -> list[float]:
    while 1:
        vect = _random_vector()
        v_len = length(vect)
        if (v_len<=0.5):
            break;
    if v_len > 0.00000000001:
        v_len = random.random()*radius / v_len
        return add(v,scale([vect[0], vect[1], vect[2]],v_len))
    else:
        return v

#------------------------------------------------------------------------------
def random_sphere(v: Sequence[float], radius: float) -> list[float]:
    while 1:
        vect = _random_vector()
        v_len = length(vect)
        if (v_len<=0.5) and (v_len!=0.0):
            break;
    return add(v,scale([vect[0], vect[1], vect[2]],2*radius/v_len))

#------------------------------------------------------------------------------
def random_vector() -> list[float]:
    while 1:
        vect = _random_vector()
        if length(vect)<=0.5:
            break;
    return scale([vect[0], vect[1], vect[2]],2.0)

#------------------------------------------------------------------------------
def add(v1: Sequence[float], v2: Sequence[float]) -> list[float]:
    return [v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]]

#------------------------------------------------------------------------------
def average(v1: Sequence[float], v2: Sequence[float]) -> list[float]:
    return [(v1[0]+v2[0])/2.0,(v1[1]+v2[1])/2.0,(v1[2]+v2[2])/2.0]

#------------------------------------------------------------------------------
def scale(v: Sequence[float], factor: float) -> list[float]:
    return [v[0]*factor,v[1]*factor,v[2]*factor]

#------------------------------------------------------------------------------
def negate(v: Sequence[float]) -> list[float]:
    return [-v[0],-v[1],-v[2]]

#------------------------------------------------------------------------------
def reverse(v: Sequence[float]) -> list[float]:
    return [ -v[0], -v[1], -v[2] ]

#------------------------------------------------------------------------------
def sub(v1: Sequence[float], v2: Sequence[float]) -> list[float]:
    return [v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]]

#------------------------------------------------------------------------------
def dot_product(v1: Sequence[float], v2: Sequence[float]) -> float:
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

#------------------------------------------------------------------------------
def cross_product(v1: Sequence[float], v2: Sequence[float]) -> list[float]:
  return [(v1[1]*v2[2]) - (v1[2]*v2[1]),
             (v1[2]*v2[0]) - (v1[0]*v2[2]),
             (v1[0]*v2[1]) - (v1[1]*v2[0])]

#------------------------------------------------------------------------------
def transform(m: Sequence[Sequence[float]], v: Sequence[float]) -> list[float]:
    return [m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
              m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
              m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]]

#------------------------------------------------------------------------------
def inverse_transform(m: Sequence[Sequence[float]], v: Sequence[float]) -> list[float]:
    return [m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2],
            m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2],
            m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2]]

#------------------------------------------------------------------------------
def multiply(m1: Sequence[Sequence[float]], m2: Sequence[Sequence[float]]) -> list[list[float]]:
   # HAVEN'T YET VERIFIED THAT THIS CONFORMS TO STANDARD DEFT
   # upd: no, it's not(fixed)

   return [[m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0],
               m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1],
               m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2]],
             [m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0],
               m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1],
               m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2]],
             [m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0],
               m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1],
               m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2]]]

#------------------------------------------------------------------------------
def get_system2(x: Sequence[float],y: Sequence[float]) -> list[list[float]]:
    z = cross_product(x,y)
    z = normalize(z)
    y = cross_product(z,x);
    y = normalize(y);
    x = normalize(x);
    return [x,y,z]

#------------------------------------------------------------------------------
def scale_system(s: Sequence[Sequence[float]], factor: float) -> list[list[float]]:
    r = []
    for a in s:
        r.append([a[0]*factor,a[1]*factor,a[2]*factor])
    return r

#------------------------------------------------------------------------------
def transpose(m: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[m[0][0], m[1][0], m[2][0]],
              [m[0][1], m[1][1], m[2][1]],
              [m[0][2], m[1][2], m[2][2]]]

#------------------------------------------------------------------------------
def transform_about_point(m: Sequence[Sequence[float]], v: Sequence[float], p: Sequence[float]) -> list[float]:
    return add(transform(m,sub(v,p)),p)

#------------------------------------------------------------------------------
def get_angle(v1: Sequence[float], v2: Sequence[float]) -> float: # v1,v2 must be unit vectors
    denom = (math.sqrt(((v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]))) *
                math.sqrt(((v2[0]*v2[0]) + (v2[1]*v2[1]) + (v2[2]*v2[2]))))
    if denom>1e-10:
        result = ( (v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]) ) / denom
    else:
        result = 0.0
    result = math.acos(result)
    return result

#------------------------------------------------------------------------------
def get_angle_formed_by(p1: Sequence[float], p2: Sequence[float], p3: Sequence[float]) -> float: # angle formed by three positions in space

    # based on code submitted by Paul Sherwood
    r1 = distance(p1,p2)
    r2 = distance(p2,p3)
    r3 = distance(p1,p3)

    small = 1.0e-10

    if (r1 + r2 - r3) < small:
        # This seems to happen occasionally for 180 angles
        theta = math.pi
    else:
        theta = math.acos( (r1*r1 + r2*r2  - r3*r3) / (2.0 * r1*r2) )
    return theta;

#------------------------------------------------------------------------------
def project(v: Sequence[float], n: Sequence[float]) -> list[float]:
    dot = v[0]*n[0] + v[1]*n[1] + v[2]*n[2]
    return [ dot * n[0], dot * n[1], dot * n[2] ]

#------------------------------------------------------------------------------
def remove_component(v: Sequence[float], n: Sequence[float]) -> list[float]:
    dot = v[0]*n[0] + v[1]*n[1] + v[2]*n[2]
    return [v[0] - dot * n[0], v[1] - dot * n[1], v[2] - dot * n[2]]

#------------------------------------------------------------------------------
def normalize(v: Sequence[float]) -> list[float]:
    vlen = math.sqrt((v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2]))
    if vlen>RSMALL4:
        return [v[0]/vlen,v[1]/vlen,v[2]/vlen]
    else:
        return get_null()

#------------------------------------------------------------------------------
def normalize_failsafe(v: Sequence[float]) -> list[float]:
    vlen = math.sqrt((v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2]))
    if vlen>RSMALL4:
        return [v[0]/vlen,v[1]/vlen,v[2]/vlen]
    else:
        return [1.0,0.0,0.0]

#------------------------------------------------------------------------------
def rotation_matrix(angle: float, axis: Sequence[float]) -> list[list[float]]:

    x=axis[0]
    y=axis[1]
    z=axis[2]

    s = math.sin(angle)
    c = math.cos(angle)

    mag = math.sqrt( x*x + y*y + z*z )

    if abs(mag)<RSMALL4:
        return get_identity()

    x = x / mag
    y = y / mag
    z = z / mag

    xx = x * x
    yy = y * y
    zz = z * z
    xy = x * y
    yz = y * z
    zx = z * x
    xs = x * s
    ys = y * s
    zs = z * s
    one_c = 1.0 - c

    return [[ (one_c * xx) + c , (one_c * xy) - zs, (one_c * zx) + ys],
              [ (one_c * xy) + zs, (one_c * yy) + c , (one_c * yz) - xs],
              [ (one_c * zx) - ys, (one_c * yz) + xs, (one_c * zz) + c ]]

#------------------------------------------------------------------------------
def transform_array(rot_mtx: Sequence[Sequence[float]], vec_array: Sequence[Sequence[float]]) -> list[list[float]]:

    '''transform_array( matrix, vector_array ) -> vector_array

    '''

    return [transform(rot_mtx, x) for x in vec_array]

#------------------------------------------------------------------------------
def translate_array(trans_vec: Sequence[float], vec_array: Sequence[Sequence[float]]) -> list[list[float]]:

    '''translate_array(trans_vec,vec_array) -> vec_array

    Adds 'mult'*'trans_vec' to each element in vec_array, and returns
    the translated vector.
    '''

    return [add(trans_vec, x) for x in vec_array]

#------------------------------------------------------------------------------
def fit_apply(fit_result: tuple[Sequence[float], Sequence[float], Sequence[Sequence[float]], float], vec_array: Sequence[Sequence[float]]) -> list[list[float]]:
    '''fit_apply(fir_result,vec_array) -> vec_array

    Applies a fit result to an array of vectors
    '''

    t1, mt2, m = fit_result[:3]
    return [add(t1, transform(m, add(mt2, x))) for x in vec_array]

#------------------------------------------------------------------------------
def fit(target_array: Sequence[Sequence[float]], source_array: Sequence[Sequence[float]]) -> tuple[list[float], list[float], list[list[float]], float]:

    '''fit(target_array, source_array) -> (t1, t2, rot_mtx, rmsd) [fit_result]

    Calculates the translation vectors and rotation matrix required
    to superimpose source_array onto target_array.  Original arrays are
    not modified.  NOTE: Currently assumes 3-dimensional coordinates

    t1,t2 are vectors from origin to centers of mass...
    '''

# Check dimensions of input arrays
    if len(target_array) != len(source_array):
        print ("Error: arrays must be of same length for RMS fitting.")
        raise ValueError
    if len(target_array[0]) != 3 or len(source_array[0]) != 3:
        print ("Error: arrays must be dimension 3 for RMS fitting.")
        raise ValueError
    nvec = len(target_array)
    ndim = 3
    maxiter = 200
    tol = 0.001

# Calculate translation vectors (center-of-mass).

    t1 = get_null()
    t2 = get_null()
    tvec1 = get_null()
    tvec2 = get_null()

    for i in range(nvec):
        for j in range(ndim):
            t1[j] = t1[j] + target_array[i][j]
            t2[j] = t2[j] + source_array[i][j]
    for j in range(ndim):
        t1[j] = t1[j] / nvec
        t2[j] = t2[j] / nvec

# Calculate correlation matrix.

    corr_mtx = []
    for i in range(ndim):
        temp_vec = []
        for j in range(ndim):
            temp_vec.append(0.0)
        corr_mtx.append(temp_vec)

    rot_mtx = []
    for i in range(ndim):
        temp_vec = []
        for j in range(ndim):
            temp_vec.append(0.0)
        rot_mtx.append(temp_vec)
    for i in range(ndim):
        rot_mtx[i][i] = 1.

    for i in range(nvec):
        for j in range(ndim):
            tvec1[j] = target_array[i][j] - t1[j]
            tvec2[j] = source_array[i][j] - t2[j]
        for j in range(ndim):
            for k in range(ndim):
                corr_mtx[j][k] = corr_mtx[j][k] + tvec2[j]*tvec1[k]

# Main iteration scheme (hardwired for 3X3 matrix, but could be extended).

    iters = 0
    while (iters < maxiter):
        iters = iters + 1
        iy = iters%ndim
        iz = (iters+1)%ndim
        sig = corr_mtx[iz][iy] - corr_mtx[iy][iz]
        gam = corr_mtx[iy][iy] + corr_mtx[iz][iz]

        sg = (sig**2 + gam**2)**0.5
        if sg != 0.0 and (abs(sig) > tol*abs(gam)):
            sg = 1.0 / sg
            for i in range(ndim):

                bb = gam*corr_mtx[iy][i] + sig*corr_mtx[iz][i]
                cc = gam*corr_mtx[iz][i] - sig*corr_mtx[iy][i]
                corr_mtx[iy][i] = bb*sg
                corr_mtx[iz][i] = cc*sg

                bb = gam*rot_mtx[iy][i] + sig*rot_mtx[iz][i]
                cc = gam*rot_mtx[iz][i] - sig*rot_mtx[iy][i]
                rot_mtx[iy][i] = bb*sg
                rot_mtx[iz][i] = cc*sg

        else:
# We have a converged rotation matrix.  Calculate RMS deviation.
            vt1 = translate_array(negate(t1),target_array)
            vt2 = translate_array(negate(t2),source_array)
            vt3 = transform_array(rot_mtx,vt2)
            rmsd = 0.0
            for i in range(nvec):
                rmsd = rmsd + distance_sq(vt1[i], vt3[i])
            rmsd = math.sqrt(rmsd/nvec)
            return(t1, t2, rot_mtx, rmsd)

# Too many iterations; something wrong.
    raise ValueError("Error: Too many iterations in RMS fit.")
