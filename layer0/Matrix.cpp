/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* Jacques Leroy (matrix inversion)
-* Thomas Malik (matrix multiplication)
-* Whoever wrote EISPACK
Z* -------------------------------------------------------------------
*/
#include"os_python.h"
#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"Vector.h"
#include"Matrix.h"
#include"MemoryDebug.h"
#include"Ortho.h"
#include"Feedback.h"
#include"Setting.h"

#include "tnt/jama_eig.h"

/* Jenarix owned types, aliases, and defines */

#define xx_os_malloc malloc
#define xx_os_free free
#define xx_os_memcpy memcpy
#define xx_os_memset memset
#define xx_fabs fabs
#define xx_sqrt sqrt
#define xx_sizeof sizeof
#define xx_float64 double
#define xx_word int
#define xx_boolean int

#ifndef XX_TRUE
#define XX_TRUE 1
#endif
#ifndef XX_FALSE
#define XX_FALSE 0
#endif
#ifndef XX_NULL
#define XX_NULL NULL
#endif

#define XX_MATRIX_STACK_STORAGE_MAX 5


/* 
   for column-major(Fortran/OpenGL) use: (col*size + row)
   for row-major(C) use: (row*size + col)
*/

#define XX_MAT(skip,row,col) (((row)*(skip)) + col)

xx_boolean xx_matrix_jacobi_solve(xx_float64 * e_vec, xx_float64 * e_val,
                                  xx_word * n_rot, const xx_float64 * input, xx_word size)


/* eigenvalues and eigenvectors for a real symmetric matrix 
   NOTE: transpose the eigenvector matrix to get eigenvectors */
{
  xx_float64 stack_A_tmp[XX_MATRIX_STACK_STORAGE_MAX * XX_MATRIX_STACK_STORAGE_MAX];
  xx_float64 stack_b_tmp[XX_MATRIX_STACK_STORAGE_MAX];
  xx_float64 stack_z_tmp[XX_MATRIX_STACK_STORAGE_MAX];
  xx_float64 *A_tmp = XX_NULL;
  xx_float64 *b_tmp = XX_NULL;
  xx_float64 *z_tmp = XX_NULL;
  xx_boolean ok = XX_TRUE;

  if(size > XX_MATRIX_STACK_STORAGE_MAX) {
    A_tmp = (xx_float64 *) xx_os_malloc(xx_sizeof(xx_float64) * size * size);
    b_tmp = (xx_float64 *) xx_os_malloc(xx_sizeof(xx_float64) * size);
    z_tmp = (xx_float64 *) xx_os_malloc(xx_sizeof(xx_float64) * size);
    if(!(A_tmp && b_tmp && z_tmp))
      ok = XX_FALSE;
  } else {
    A_tmp = stack_A_tmp;
    b_tmp = stack_b_tmp;
    z_tmp = stack_z_tmp;
  }

  if(ok) {

    xx_os_memset(e_vec, 0, xx_sizeof(xx_float64) * size * size);
    xx_os_memcpy(A_tmp, input, xx_sizeof(xx_float64) * size * size);

    {
      xx_word p;
      for(p = 0; p < size; p++) {
        e_vec[XX_MAT(size, p, p)] = 1.0;
        e_val[p] = b_tmp[p] = A_tmp[XX_MAT(size, p, p)];
        z_tmp[p] = 0.0;
      }
      *n_rot = 0;
    }
    {
      xx_word i;
      for(i = 0; i < 50; i++) {
        xx_float64 thresh;

        {
          xx_word p, q;
          xx_float64 sm = 0.0;
          for(p = 0; p < (size - 1); p++) {
            for(q = p + 1; q < size; q++) {
              sm += xx_fabs(A_tmp[XX_MAT(size, p, q)]);
            }
          }
          if(sm == 0.0)
            break;
          if(i < 3) {
            thresh = 0.2 * sm / (size * size);
          } else {
            thresh = 0.0;
          }
        }

        {
          xx_word p, q;
          xx_float64 g;

          for(p = 0; p < (size - 1); p++) {
            for(q = p + 1; q < size; q++) {
              g = 100.0 * xx_fabs(A_tmp[XX_MAT(size, p, q)]);
              if((i > 3) &&
                 ((xx_fabs(e_val[p]) + g) == xx_fabs(e_val[p])) &&
                 ((xx_fabs(e_val[q]) + g) == xx_fabs(e_val[q]))) {
                A_tmp[XX_MAT(size, p, q)] = 0.0;
              } else if(xx_fabs(A_tmp[XX_MAT(size, p, q)]) > thresh) {
                xx_float64 t;
                xx_float64 h = e_val[q] - e_val[p];
                if((xx_fabs(h) + g) == xx_fabs(h)) {
                  t = A_tmp[XX_MAT(size, p, q)] / h;
                } else {
                  xx_float64 theta = 0.5 * h / A_tmp[XX_MAT(size, p, q)];
                  t = 1.0 / (xx_fabs(theta) + xx_sqrt(1.0 + theta * theta));
                  if(theta < 0.0)
                    t = -t;
                }
                {
                  xx_float64 c = 1.0 / xx_sqrt(1.0 + t * t);
                  xx_float64 s = t * c;
                  xx_float64 tau = s / (1.0 + c);

                  h = t * A_tmp[XX_MAT(size, p, q)];
                  z_tmp[p] -= h;
                  z_tmp[q] += h;
                  e_val[p] -= h;
                  e_val[q] += h;
                  A_tmp[XX_MAT(size, p, q)] = 0.0;
                  {
                    xx_word j;
                    for(j = 0; j < p; j++) {
                      g = A_tmp[XX_MAT(size, j, p)];
                      h = A_tmp[XX_MAT(size, j, q)];
                      A_tmp[XX_MAT(size, j, p)] = g - s * (h + g * tau);
                      A_tmp[XX_MAT(size, j, q)] = h + s * (g - h * tau);
                    }
                    for(j = p + 1; j < q; j++) {
                      g = A_tmp[XX_MAT(size, p, j)];
                      h = A_tmp[XX_MAT(size, j, q)];
                      A_tmp[XX_MAT(size, p, j)] = g - s * (h + g * tau);
                      A_tmp[XX_MAT(size, j, q)] = h + s * (g - h * tau);
                    }
                    for(j = q + 1; j < size; j++) {
                      g = A_tmp[XX_MAT(size, p, j)];
                      h = A_tmp[XX_MAT(size, q, j)];
                      A_tmp[XX_MAT(size, p, j)] = g - s * (h + g * tau);
                      A_tmp[XX_MAT(size, q, j)] = h + s * (g - h * tau);
                    }
                    for(j = 0; j < size; j++) {
                      g = e_vec[XX_MAT(size, j, p)];
                      h = e_vec[XX_MAT(size, j, q)];
                      e_vec[XX_MAT(size, j, p)] = g - s * (h + g * tau);
                      e_vec[XX_MAT(size, j, q)] = h + s * (g - h * tau);
                    }
                  }
                  (*n_rot)++;
                }
              }
            }
          }
          for(p = 0; p < size; p++) {
            b_tmp[p] += z_tmp[p];
            e_val[p] = b_tmp[p];
            z_tmp[p] = 0.0;
          }
        }
      }
    }
  }
  if(A_tmp && (A_tmp != stack_A_tmp))
    xx_os_free(A_tmp);
  if(b_tmp && (b_tmp != stack_b_tmp))
    xx_os_free(b_tmp);
  if(z_tmp && (z_tmp != stack_z_tmp))
    xx_os_free(z_tmp);
  return ok;
}

static xx_word xx_matrix_decompose(xx_float64 * matrix, xx_word size,
                                   xx_word * permute, xx_word * parity)
{

  xx_boolean ok = XX_TRUE;
  xx_float64 stack_storage[XX_MATRIX_STACK_STORAGE_MAX];
  xx_float64 *storage = XX_NULL;

  if(size > XX_MATRIX_STACK_STORAGE_MAX) {
    storage = (xx_float64 *) xx_os_malloc(xx_sizeof(xx_float64) * size);
    if(!storage)
      ok = false;
  } else {
    storage = stack_storage;
  }

  *parity = 1;

  if(ok) {
    xx_word i, j;
    for(i = 0; i < size; i++) {
      xx_float64 max_abs_value = 0.0;
      for(j = 0; j < size; j++) {
        xx_float64 test_abs_value = xx_fabs(matrix[XX_MAT(size, i, j)]);
        if(max_abs_value < test_abs_value)
          max_abs_value = test_abs_value;
      }
      if(max_abs_value == 0.0) {
        ok = XX_FALSE;          /* singular matrix -- no inverse */
        break;
      }
      storage[i] = 1.0 / max_abs_value;
    }
  }
  if(ok) {
    xx_word i, j, k, i_max = 0;

    for(j = 0; j < size; j++) {

      for(i = 0; i < j; i++) {
        xx_float64 sum = matrix[XX_MAT(size, i, j)];
        for(k = 0; k < i; k++) {
          sum = sum - matrix[XX_MAT(size, i, k)] * matrix[XX_MAT(size, k, j)];
        }
        matrix[XX_MAT(size, i, j)] = sum;
      }

      {
        xx_float64 max_product = 0.0;
        for(i = j; i < size; i++) {
          xx_float64 sum = matrix[XX_MAT(size, i, j)];
          for(k = 0; k < j; k++) {
            sum = sum - matrix[XX_MAT(size, i, k)] * matrix[XX_MAT(size, k, j)];
          }
          matrix[XX_MAT(size, i, j)] = sum;
          {
            xx_float64 test_product = storage[i] * xx_fabs(sum);
            if(max_product <= test_product) {
              max_product = test_product;
              i_max = i;
            }
          }
        }
      }

      if(j != i_max) {
        for(k = 0; k < size; k++) {
          xx_float64 tmp = matrix[XX_MAT(size, i_max, k)];
          matrix[XX_MAT(size, i_max, k)] = matrix[XX_MAT(size, j, k)];
          matrix[XX_MAT(size, j, k)] = tmp;
        }
        *parity = (0 - *parity);
        storage[i_max] = storage[j];
      }

      permute[j] = i_max;
      if(matrix[XX_MAT(size, j, j)] == 0.0) {
        /* here we have a choice: */
        /* or (B): fail outright */
        ok = false;
        break;
      }

      if(j != (size - 1)) {
        xx_float64 tmp = 1.0 / matrix[XX_MAT(size, j, j)];
        for(i = j + 1; i < size; i++) {
          matrix[XX_MAT(size, i, j)] *= tmp;
        }
      }
    }
  }
  if(storage && (storage != stack_storage))
    xx_os_free(storage);
  return ok;
}

static void xx_matrix_back_substitute(xx_float64 * result, xx_float64 * decomp,
                                      xx_word size, xx_word * permute)
{
  {
    xx_word i, ii;
    ii = -1;
    for(i = 0; i < size; i++) {
      xx_word p = permute[i];
      xx_float64 sum = result[p];
      result[p] = result[i];
      if(ii >= 0) {
        xx_word j;
        for(j = ii; j < i; j++) {
          sum = sum - decomp[XX_MAT(size, i, j)] * result[j];
        }
      } else if(sum != 0.0) {
        ii = i;
      }
      result[i] = sum;
    }
  }

  {
    xx_word i, j;
    for(i = size - 1; i >= 0; i--) {
      xx_float64 sum = result[i];
      for(j = i + 1; j < size; j++) {
        sum = sum - decomp[XX_MAT(size, i, j)] * result[j];
      }
      result[i] = sum / decomp[XX_MAT(size, i, i)];
    }
  }
}

xx_boolean xx_matrix_invert(xx_float64 * result, const xx_float64 * input, xx_word size)
{
  xx_float64 stack_mat_tmp[XX_MATRIX_STACK_STORAGE_MAX * XX_MATRIX_STACK_STORAGE_MAX];
  xx_float64 stack_dbl_tmp[XX_MATRIX_STACK_STORAGE_MAX];
  xx_word stack_int_tmp[XX_MATRIX_STACK_STORAGE_MAX];
  xx_float64 *mat_tmp = XX_NULL;
  xx_float64 *dbl_tmp = XX_NULL;
  xx_word *int_tmp = XX_NULL;
  xx_word parity = 0;
  xx_word ok = XX_TRUE;

  if(size > XX_MATRIX_STACK_STORAGE_MAX) {
    mat_tmp = (xx_float64 *) xx_os_malloc(xx_sizeof(xx_float64) * size * size);
    dbl_tmp = (xx_float64 *) xx_os_malloc(xx_sizeof(xx_float64) * size);
    int_tmp = (xx_word *) xx_os_malloc(xx_sizeof(xx_word) * size);
    if(!(mat_tmp && dbl_tmp && int_tmp))
      ok = XX_FALSE;
  } else {
    mat_tmp = stack_mat_tmp;
    dbl_tmp = stack_dbl_tmp;
    int_tmp = stack_int_tmp;
  }

  if(ok) {
    ok = XX_FALSE;
    memcpy(mat_tmp, input, xx_sizeof(xx_float64) * size * size);

    if(xx_matrix_decompose(mat_tmp, size, int_tmp, &parity)) {
      xx_word i, j;
      for(j = 0; j < size; j++) {
        memset(dbl_tmp, 0, xx_sizeof(xx_float64) * size);
        dbl_tmp[j] = 1.0;
        xx_matrix_back_substitute(dbl_tmp, mat_tmp, size, int_tmp);
        for(i = 0; i < size; i++) {
          result[XX_MAT(size, i, j)] = dbl_tmp[i];
        }
      }
      ok = XX_TRUE;
    }
  }
  /* clean up */
  if(mat_tmp && (mat_tmp != stack_mat_tmp))
    xx_os_free(mat_tmp);
  if(dbl_tmp && (dbl_tmp != stack_dbl_tmp))
    xx_os_free(dbl_tmp);
  if(int_tmp && (int_tmp != stack_int_tmp))
    xx_os_free(int_tmp);
  return ok;
}

#undef XX_MAT

int MatrixInvTransformExtentsR44d3f(const double *matrix,
                                    const float *old_min, const float *old_max,
                                    float *new_min, float *new_max)
{
  /* just brute-forcing this for now... */
  int a;
  int c;

  double inp_min[3], inp_max[3];
  double out_min[3], out_max[3];
  double inp_tst[3], out_tst[3];

  if(!matrix)
    return 0;

  copy3f3d(old_min, inp_min);
  copy3f3d(old_max, inp_max);

  for(c = 0; c < 8; c++) {
    inp_tst[0] = c & 0x1 ? inp_min[0] : inp_max[0];
    inp_tst[1] = c & 0x2 ? inp_min[1] : inp_max[1];
    inp_tst[2] = c & 0x4 ? inp_min[2] : inp_max[2];

    inverse_transform44d3d(matrix, inp_tst, out_tst);
    if(!c) {
      copy3d(out_tst, out_max);
      copy3d(out_tst, out_min);
    } else {
      for(a = 0; a < 3; a++) {
        if(out_min[a] > out_tst[a])
          out_min[a] = out_tst[a];
        if(out_max[a] < out_tst[a])
          out_max[a] = out_tst[a];
      }
    }
  }
  copy3d3f(out_min, new_min);
  copy3d3f(out_max, new_max);
  return 1;
}

int MatrixTransformExtentsR44d3f(const double *matrix,
                                 const float *old_min, const float *old_max,
                                 float *new_min, float *new_max)
{
  /* just brute-forcing this for now... */
  int a;
  int c;

  double inp_min[3], inp_max[3];
  double out_min[3], out_max[3];
  double inp_tst[3], out_tst[3];

  if(!matrix)
    return 0;

  copy3f3d(old_min, inp_min);
  copy3f3d(old_max, inp_max);

  for(c = 0; c < 8; c++) {
    inp_tst[0] = c & 0x1 ? inp_min[0] : inp_max[0];
    inp_tst[1] = c & 0x2 ? inp_min[1] : inp_max[1];
    inp_tst[2] = c & 0x4 ? inp_min[2] : inp_max[2];

    transform44d3d(matrix, inp_tst, out_tst);
    if(!c) {
      copy3d(out_tst, out_max);
      copy3d(out_tst, out_min);
    } else {
      for(a = 0; a < 3; a++) {
        if(out_min[a] > out_tst[a])
          out_min[a] = out_tst[a];
        if(out_max[a] < out_tst[a])
          out_max[a] = out_tst[a];
      }
    }
  }
  copy3d3f(out_min, new_min);
  copy3d3f(out_max, new_max);
  return 1;
}

int MatrixInvertC44f(const float *m, float *out)
{

  /* This routine included in PyMOL under the terms of the 
   * MIT consortium license for Brian Paul's Mesa, from which it was derived. */

  /* MESA comments:
   * Compute inverse of 4x4 transformation matrix.
   * Code contributed by Jacques Leroy jle@star.be
   * Return GL_TRUE for success, GL_FALSE for failure (singular matrix)
   */

  /* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { float *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(m,r,c) (m)[(c)*4+(r)]

  float wtmp[4][8];
  float m0, m1, m2, m3, s;
  float *r0, *r1, *r2, *r3;

  r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

  r0[0] = MAT(m, 0, 0), r0[1] = MAT(m, 0, 1),
    r0[2] = MAT(m, 0, 2), r0[3] = MAT(m, 0, 3),
    r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0F,
    r1[0] = MAT(m, 1, 0), r1[1] = MAT(m, 1, 1),
    r1[2] = MAT(m, 1, 2), r1[3] = MAT(m, 1, 3),
    r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0F,
    r2[0] = MAT(m, 2, 0), r2[1] = MAT(m, 2, 1),
    r2[2] = MAT(m, 2, 2), r2[3] = MAT(m, 2, 3),
    r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0F,
    r3[0] = MAT(m, 3, 0), r3[1] = MAT(m, 3, 1),
    r3[2] = MAT(m, 3, 2), r3[3] = MAT(m, 3, 3), r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0F;

  /* choose pivot - or die */
  if(fabs(r3[0]) > fabs(r2[0]))
    SWAP_ROWS(r3, r2);
  if(fabs(r2[0]) > fabs(r1[0]))
    SWAP_ROWS(r2, r1);
  if(fabs(r1[0]) > fabs(r0[0]))
    SWAP_ROWS(r1, r0);
  if(0.0F == r0[0])
    return 0;

  /* eliminate first variable     */
  m1 = r1[0] / r0[0];
  m2 = r2[0] / r0[0];
  m3 = r3[0] / r0[0];
  s = r0[1];
  r1[1] -= m1 * s;
  r2[1] -= m2 * s;
  r3[1] -= m3 * s;
  s = r0[2];
  r1[2] -= m1 * s;
  r2[2] -= m2 * s;
  r3[2] -= m3 * s;
  s = r0[3];
  r1[3] -= m1 * s;
  r2[3] -= m2 * s;
  r3[3] -= m3 * s;
  s = r0[4];
  if(s != 0.0F) {
    r1[4] -= m1 * s;
    r2[4] -= m2 * s;
    r3[4] -= m3 * s;
  }
  s = r0[5];
  if(s != 0.0F) {
    r1[5] -= m1 * s;
    r2[5] -= m2 * s;
    r3[5] -= m3 * s;
  }
  s = r0[6];
  if(s != 0.0F) {
    r1[6] -= m1 * s;
    r2[6] -= m2 * s;
    r3[6] -= m3 * s;
  }
  s = r0[7];
  if(s != 0.0F) {
    r1[7] -= m1 * s;
    r2[7] -= m2 * s;
    r3[7] -= m3 * s;
  }

  /* choose pivot - or die */
  if(fabs(r3[1]) > fabs(r2[1]))
    SWAP_ROWS(r3, r2);
  if(fabs(r2[1]) > fabs(r1[1]))
    SWAP_ROWS(r2, r1);
  if(0.0F == r1[1])
    return 0;

  /* eliminate second variable */
  m2 = r2[1] / r1[1];
  m3 = r3[1] / r1[1];
  r2[2] -= m2 * r1[2];
  r3[2] -= m3 * r1[2];
  r2[3] -= m2 * r1[3];
  r3[3] -= m3 * r1[3];
  s = r1[4];
  if(0.0F != s) {
    r2[4] -= m2 * s;
    r3[4] -= m3 * s;
  }
  s = r1[5];
  if(0.0F != s) {
    r2[5] -= m2 * s;
    r3[5] -= m3 * s;
  }
  s = r1[6];
  if(0.0F != s) {
    r2[6] -= m2 * s;
    r3[6] -= m3 * s;
  }
  s = r1[7];
  if(0.0F != s) {
    r2[7] -= m2 * s;
    r3[7] -= m3 * s;
  }

  /* choose pivot - or die */
  if(fabs(r3[2]) > fabs(r2[2]))
    SWAP_ROWS(r3, r2);
  if(0.0F == r2[2])
    return 0;

  /* eliminate third variable */
  m3 = r3[2] / r2[2];
  r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
    r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

  /* last check */
  if(0.0F == r3[3])
    return 0;

  s = 1.0F / r3[3];             /* now back substitute row 3 */
  r3[4] *= s;
  r3[5] *= s;
  r3[6] *= s;
  r3[7] *= s;

  m2 = r2[3];                   /* now back substitute row 2 */
  s = 1.0F / r2[2];
  r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
    r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
  m1 = r1[3];
  r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1, r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
  m0 = r0[3];
  r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0, r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

  m1 = r1[2];                   /* now back substitute row 1 */
  s = 1.0F / r1[1];
  r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
    r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
  m0 = r0[2];
  r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0, r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

  m0 = r0[1];                   /* now back substitute row 0 */
  s = 1.0F / r0[0];
  r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
    r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

  MAT(out, 0, 0) = r0[4];
  MAT(out, 0, 1) = r0[5], MAT(out, 0, 2) = r0[6];
  MAT(out, 0, 3) = r0[7], MAT(out, 1, 0) = r1[4];
  MAT(out, 1, 1) = r1[5], MAT(out, 1, 2) = r1[6];
  MAT(out, 1, 3) = r1[7], MAT(out, 2, 0) = r2[4];
  MAT(out, 2, 1) = r2[5], MAT(out, 2, 2) = r2[6];
  MAT(out, 2, 3) = r2[7], MAT(out, 3, 0) = r3[4];
  MAT(out, 3, 1) = r3[5], MAT(out, 3, 2) = r3[6];
  MAT(out, 3, 3) = r3[7];

  return 1;

#undef MAT
#undef SWAP_ROWS
}

/*========================================================================*/
void MatrixTransformTTTfN3f(unsigned int n, float *q, const float *m, const float *p)
{
  const float m0 = m[0], m4 = m[4], m8 = m[8], m12 = m[12];
  const float m1 = m[1], m5 = m[5], m9 = m[9], m13 = m[13];
  const float m2 = m[2], m6 = m[6], m10 = m[10], m14 = m[14];
  const float m3 = m[3], m7 = m[7], m11 = m[11];
  float p0, p1, p2;
  while(n--) {
    p0 = *(p++) + m12;
    p1 = *(p++) + m13;
    p2 = *(p++) + m14;
    *(q++) = m0 * p0 + m1 * p1 + m2 * p2 + m3;
    *(q++) = m4 * p0 + m5 * p1 + m6 * p2 + m7;
    *(q++) = m8 * p0 + m9 * p1 + m10 * p2 + m11;
  }
}


/*========================================================================*/
void MatrixTransformR44fN3f(unsigned int n, float *q, const float *m, const float *p)
{
  const float m0 = m[0], m4 = m[4], m8 = m[8];
  const float m1 = m[1], m5 = m[5], m9 = m[9];
  const float m2 = m[2], m6 = m[6], m10 = m[10];
  const float m3 = m[3], m7 = m[7], m11 = m[11];
  float p0, p1, p2;
  while(n--) {
    p0 = *(p++);
    p1 = *(p++);
    p2 = *(p++);
    *(q++) = m0 * p0 + m1 * p1 + m2 * p2 + m3;
    *(q++) = m4 * p0 + m5 * p1 + m6 * p2 + m7;
    *(q++) = m8 * p0 + m9 * p1 + m10 * p2 + m11;
  }
}


/*========================================================================*/
void MatrixGetRotationC44f(float *m44, float angle,
                           float x, float y, float z)
{
  float m33[9];
  rotation_matrix3f(angle, x, y, z, m33);
  m44[0] = m33[0];
  m44[1] = m33[3];
  m44[2] = m33[6];
  m44[3] = 0.0F;
  m44[4] = m33[1];
  m44[5] = m33[4];
  m44[6] = m33[7];
  m44[7] = 0.0F;
  m44[8] = m33[2];
  m44[9] = m33[5];
  m44[10] = m33[8];
  m44[11] = 0.0F;
  m44[12] = 0.0F;
  m44[13] = 0.0F;
  m44[14] = 0.0F;
  m44[15] = 1.0F;
}


/*========================================================================*/
void MatrixRotateC44f(float *m, float angle, float x, float y, float z)
{
  float m33[9];
  float m44[16];
  rotation_matrix3f(angle, x, y, z, m33);
  m44[0] = m33[0];
  m44[1] = m33[1];
  m44[2] = m33[2];
  m44[3] = 0.0F;
  m44[4] = m33[3];
  m44[5] = m33[4];
  m44[6] = m33[5];
  m44[7] = 0.0F;
  m44[8] = m33[6];
  m44[9] = m33[7];
  m44[10] = m33[8];
  m44[11] = 0.0F;
  m44[12] = 0.0F;
  m44[13] = 0.0F;
  m44[14] = 0.0F;
  m44[15] = 1.0F;
  MatrixMultiplyC44f(m44, m);
}


/*========================================================================*/
void MatrixTranslateC44f(float *m, float x, float y, float z)
{
  m[12] = m[0] * x + m[4] * y + m[8] * z + m[12];
  m[13] = m[1] * x + m[5] * y + m[9] * z + m[13];
  m[14] = m[2] * x + m[6] * y + m[10] * z + m[14];
  m[15] = m[3] * x + m[7] * y + m[11] * z + m[15];
}


/*========================================================================*/
void MatrixMultiplyC44f(const float *b, float *m)
{
  /* This routine included in PyMOL under the terms of the 
   * MIT consortium license for Brian Paul's Mesa, from which it was derived. */
  /* original author: Thomas Malik */

  int i;

#define A(row,col)  m[(col<<2)+row]
#define B(row,col)  b[(col<<2)+row]
#define P(row,col)  m[(col<<2)+row]

  /* i-te Zeile */
  for(i = 0; i < 4; i++) {
    float ai0 = A(i, 0), ai1 = A(i, 1), ai2 = A(i, 2), ai3 = A(i, 3);
    P(i, 0) = ai0 * B(0, 0) + ai1 * B(1, 0) + ai2 * B(2, 0) + ai3 * B(3, 0);
    P(i, 1) = ai0 * B(0, 1) + ai1 * B(1, 1) + ai2 * B(2, 1) + ai3 * B(3, 1);
    P(i, 2) = ai0 * B(0, 2) + ai1 * B(1, 2) + ai2 * B(2, 2) + ai3 * B(3, 2);
    P(i, 3) = ai0 * B(0, 3) + ai1 * B(1, 3) + ai2 * B(2, 3) + ai3 * B(3, 3);
  }

#undef A
#undef B
#undef P
}


/*========================================================================*/
void MatrixTransformC44f3f(const float *m, const float *q, float *p)
{
  float q0 = *q, q1 = *(q + 1), q2 = *(q + 2);
  p[0] = m[0] * q0 + m[4] * q1 + m[8] * q2 + m[12];
  p[1] = m[1] * q0 + m[5] * q1 + m[9] * q2 + m[13];
  p[2] = m[2] * q0 + m[6] * q1 + m[10] * q2 + m[14];
}


/*========================================================================*/
void MatrixTransformC44f4f(const float *m, const float *q, float *p)
{
  float q0 = *q, q1 = *(q + 1), q2 = *(q + 2);
  p[0] = m[0] * q0 + m[4] * q1 + m[8] * q2 + m[12];
  p[1] = m[1] * q0 + m[5] * q1 + m[9] * q2 + m[13];
  p[2] = m[2] * q0 + m[6] * q1 + m[10] * q2 + m[14];
  p[3] = m[3] * q0 + m[7] * q1 + m[11] * q2 + m[15];
}


/*========================================================================*/
// same as transform44f3fas33f3f
void MatrixInvTransformC44fAs33f3f(const float *m, const float *q, float *p)
{
  /* multiplying a column major rotation matrix as row-major will
   * give the inverse rotation */
  float q0 = *q, q1 = *(q + 1), q2 = *(q + 2);
  p[0] = m[0] * q0 + m[1] * q1 + m[2] * q2;
  p[1] = m[4] * q0 + m[5] * q1 + m[6] * q2;
  p[2] = m[8] * q0 + m[9] * q1 + m[10] * q2;
}


/*========================================================================*/
void MatrixTransformC44fAs33f3f(const float *m, const float *q, float *p)
{
  float q0 = *q, q1 = *(q + 1), q2 = *(q + 2);
  p[0] = m[0] * q0 + m[4] * q1 + m[8] * q2;
  p[1] = m[1] * q0 + m[5] * q1 + m[9] * q2;
  p[2] = m[2] * q0 + m[6] * q1 + m[10] * q2;
}


/*========================================================================*/
float MatrixGetRMS(PyMOLGlobals * G, int n, const float *v1, const float *v2, float *wt)
{
  /* Just Compute RMS given current coordinates */

  const float *vv1, *vv2;
  float err, etmp, tmp;
  int a, c;
  float sumwt = 0.0F;

  if(wt) {
    for(c = 0; c < n; c++)
      if(wt[c] != 0.0F) {
        sumwt = sumwt + wt[c];
      }
  } else {
    for(c = 0; c < n; c++)
      sumwt += 1.0F;
  }
  err = 0.0F;
  vv1 = v1;
  vv2 = v2;
  for(c = 0; c < n; c++) {
    etmp = 0.0F;
    for(a = 0; a < 3; a++) {
      tmp = (vv2[a] - vv1[a]);
      etmp += tmp * tmp;
    }
    if(wt)
      err += wt[c] * etmp;
    else
      err += etmp;
    vv1 += 3;
    vv2 += 3;
  }
  err = err / sumwt;
  err = (float) sqrt1f(err);

  if(fabs(err) < R_SMALL4)
    err = 0.0F;

  return (err);
}


/*========================================================================*/
float MatrixFitRMSTTTf(PyMOLGlobals * G, int n, const float *v1, const float *v2, const float *wt,
                       float *ttt)
{
  /*
     Subroutine to do the actual RMS fitting of two sets of vector coordinates
     This routine does not rotate the actual coordinates, but instead returns 
     the RMS fitting value, along with the center-of-mass translation vectors 
     T1 and T2 and the rotation matrix M, which rotates the translated 
     coordinates of molecule 2 onto the translated coordinates of molecule 1.
   */

  double m[3][3], aa[3][3];
  double sumwt, tol;
  int a, b, c, maxiter;
  double t1[3], t2[3];

  /* special case: just one atom pair */

  if(n == 1) {
    if(ttt) {
      for(a = 1; a < 11; a++)
        ttt[a] = 0.0F;
      for(a = 0; a < 12; a += 5)
        ttt[a] = 1.0F;
      for(a = 0; a < 3; a++)
        ttt[a + 12] = v2[a] - v1[a];
    }
    return 0.0F;
  }

  /* Initialize arrays. */

  for(a = 0; a < 3; a++) {
    for(b = 0; b < 3; b++) {
      aa[a][b] = 0.0;
    }
    t1[a] = 0.0;
    t2[a] = 0.0;
  }

  sumwt = 0.0F;
  tol = SettingGetGlobal_f(G, cSetting_fit_tolerance);
  maxiter = SettingGetGlobal_i(G, cSetting_fit_iterations);

  /* Calculate center-of-mass vectors */

  {
    const float *vv1 = v1, *vv2 = v2;

    if(wt) {
      for(c = 0; c < n; c++) {
        for(a = 0; a < 3; a++) {
          t1[a] += wt[c] * vv1[a];
          t2[a] += wt[c] * vv2[a];
        }
        if(wt[c] != 0.0F) {
          sumwt = sumwt + wt[c];
        } else {
          sumwt = sumwt + 1.0F; /* WHAT IS THIS? */
        }
        vv1 += 3;
        vv2 += 3;
      }
    } else {
      for(c = 0; c < n; c++) {
        for(a = 0; a < 3; a++) {
          t1[a] += vv1[a];
          t2[a] += vv2[a];
        }
        sumwt += 1.0F;
        vv1 += 3;
        vv2 += 3;
      }
    }
    if(sumwt == 0.0F)
      sumwt = 1.0F;
    for(a = 0; a < 3; a++) {
      t1[a] /= sumwt;
      t2[a] /= sumwt;
    }
  }

  {
    /* Calculate correlation matrix */
    double x[3], xx[3];
    const float *vv1 = v1, *vv2 = v2;
    for(c = 0; c < n; c++) {
      if(wt) {
        for(a = 0; a < 3; a++) {
          x[a] = wt[c] * (vv1[a] - t1[a]);
          xx[a] = wt[c] * (vv2[a] - t2[a]);
        }
      } else {
        for(a = 0; a < 3; a++) {
          x[a] = vv1[a] - t1[a];
          xx[a] = vv2[a] - t2[a];
        }
      }
      for(a = 0; a < 3; a++)
        for(b = 0; b < 3; b++)
          aa[a][b] = aa[a][b] + xx[a] * x[b];
      vv1 += 3;
      vv2 += 3;
    }
  }

  if(n > 1) {
    int got_kabsch = false;
    int fit_kabsch = SettingGetGlobal_i(G, cSetting_fit_kabsch);
    if(fit_kabsch) {

      /* WARNING: Kabsch isn't numerically stable */

      /* Kabsch as per

         http://en.wikipedia.org/wiki/Kabsch_algorithm

         minimal RMS matrix is (AtA)^(1/2) * A_inverse, where

         Aij =  Pki Qkj

         assuming Pki and Qkj are centered about the same origin.

       */

      /* NOTE: This Kabsch implementation only works with 4 or more atoms */

      double At[3][3], AtA[3][3];

      /* compute At and At * A */

      transpose33d33d((double *) (void *) aa, (double *) (void *) At);

      multiply33d33d((double *) (void *) At, (double *) (void *) aa,
                     (double *) (void *) AtA);

      /* solve A*At (a real symmetric matrix) */

      {
        double e_vec[3][3], e_val[3];
        xx_word n_rot;

        if(xx_matrix_jacobi_solve((xx_float64 *) (void *) e_vec,
                                  (xx_float64 *) (void *) e_val,
                                  &n_rot, (xx_float64 *) (void *) AtA, 3)) {
          double V[3][3], Vt[3][3];

          /* Kabsch requires non-negative eigenvalues (since sqrt is taken) */

          if((e_val[0] >= 0.0) && (e_val[1] >= 0.0) && (e_val[2] >= 0.0)) {
            switch (fit_kabsch) {
            case 1:
              {

                /* original Kabsch performs an unnecessary matrix inversion */

                /* rot matrix = (AtA)^(1/2) * A^(-1) */

                double rootD[3][3], sqrtAtA[3][3], Ai[3][3];

                for(a = 0; a < 3; a++)
                  for(b = 0; b < 3; b++) {
                    if(a == b) {
                      rootD[a][b] = sqrt1d(e_val[a]);
                    } else {
                      rootD[a][b] = 0.0;
                    }
                    V[a][b] = e_vec[a][b];
                    Vt[a][b] = e_vec[b][a];
                  }

                multiply33d33d((double *) (void *) rootD, (double *) (void *) Vt,
                               (double *) (void *) sqrtAtA);
                multiply33d33d((double *) (void *) V, (double *) (void *) sqrtAtA,
                               (double *) (void *) sqrtAtA);
                /* compute Ai */

                if(xx_matrix_invert((xx_float64 *) (void *) Ai,
                                    (xx_float64 *) (void *) aa, 3)) {

                  /* now compute the rotation matrix  = (AtA)^(1/2) * Ai */

                  multiply33d33d((double *) (void *) sqrtAtA,
                                 (double *) (void *) Ai, (double *) (void *) m);

                  got_kabsch = true;

                  {             /* is the rotation matrix left-handed? Then swap the
                                   recomposition so as to avoid the reflection */
                    double cp[3], dp;
                    cross_product3d((double *) (void *) m, (double *) (void *) m + 3, cp);
                    dp = dot_product3d((double *) (void *) m + 6, cp);
                    if((1.0 - fabs(dp)) > R_SMALL4) {   /* not orthonormal? */
                      got_kabsch = false;

                      PRINTFB(G, FB_Matrix, FB_Warnings)
                        "Matrix-Warning: Kabsch matrix not orthonormal: falling back on iteration.\n"
                        ENDFB(G);

                    } else if(dp < 0.0F) {
                      multiply33d33d((double *) (void *) rootD, (double *) (void *) V,
                                     (double *) (void *) sqrtAtA);
                      multiply33d33d((double *) (void *) Vt, (double *) (void *) sqrtAtA,
                                     (double *) (void *) sqrtAtA);
                      multiply33d33d((double *) (void *) sqrtAtA, (double *) (void *) Ai,
                                     (double *) (void *) m);
                    }
                  }
                } else {
                  PRINTFB(G, FB_Matrix, FB_Warnings)
                    "Matrix-Warning: Kabsch matrix inversion failed: falling back on iteration.\n"
                    ENDFB(G);
                }
              }
              break;
            case 2:
              {
                /* improved Kabsch skips the matrix inversion */

                if((e_val[0] > R_SMALL8) &&
                   (e_val[1] > R_SMALL8) && (e_val[2] > R_SMALL8)) {

                  /* inv rot matrix = U Vt = A V D^(-1/2) Vt
                   * where U from SVD of A = U S Vt 
                   * and where U is known to be = A V D^(-1/2)
                   * where D are eigenvalues of AtA */

                  double invRootD[3][3], Mi[3][3];

                  for(a = 0; a < 3; a++)
                    for(b = 0; b < 3; b++) {
                      if(a == b) {
                        invRootD[a][b] = 1.0 / sqrt1d(e_val[a]);
                      } else {
                        invRootD[a][b] = 0.0;
                      }
                      V[a][b] = e_vec[a][b];
                      Vt[a][b] = e_vec[b][a];
                    }

                  /* now compute the rotation matrix directly  */

                  multiply33d33d((double *) (void *) invRootD, (double *) (void *) Vt,
                                 (double *) (void *) Mi);
                  multiply33d33d((double *) (void *) V, (double *) (void *) Mi,
                                 (double *) (void *) Mi);
                  multiply33d33d((double *) (void *) aa, (double *) (void *) Mi,
                                 (double *) (void *) Mi);

                  got_kabsch = true;

                  {             /* cis the rotation matrix left-handed? Then swap the
                                   recomposition so as to avoid the reflection */
                    double cp[3], dp;
                    cross_product3d((double *) (void *) Mi, (double *) (void *) Mi + 3,
                                    cp);
                    dp = dot_product3d((double *) (void *) Mi + 6, cp);
                    if((1.0 - fabs(dp)) > R_SMALL4) {   /* not orthonormal? */
                      got_kabsch = false;

                      PRINTFB(G, FB_Matrix, FB_Warnings)
                        "Matrix-Warning: Kabsch matrix not orthonormal: falling back on iteration.\n"
                        ENDFB(G);

                    } else if(dp < 0.0F) {
                      multiply33d33d((double *) (void *) invRootD, (double *) (void *) V,
                                     (double *) (void *) Mi);
                      multiply33d33d((double *) (void *) Vt, (double *) (void *) Mi,
                                     (double *) (void *) Mi);
                      multiply33d33d((double *) (void *) aa, (double *) (void *) Mi,
                                     (double *) (void *) Mi);
                    }
                  }

                  if(got_kabsch) {      /* transpose to get the inverse */
                    transpose33d33d((double *) (void *) Mi, (double *) (void *) m);
                  }

                } else {
                  PRINTFB(G, FB_Matrix, FB_Warnings)
                    "Matrix-Warning: Kabsch matrix degenerate: falling back on iteration.\n"
                    ENDFB(G);
                }
              }
              break;
            }

            /* validate the result */

            if(got_kabsch) {

              if((fabs(length3d(m[0]) - 1.0) < R_SMALL4) &&
                 (fabs(length3d(m[1]) - 1.0) < R_SMALL4) &&
                 (fabs(length3d(m[2]) - 1.0) < R_SMALL4)) {

                recondition33d((double *) (void *) m);

              } else {
                got_kabsch = false;

                PRINTFB(G, FB_Matrix, FB_Warnings)
                  "Matrix-Warning: Kabsch matrix not normal: falling back on iteration.\n"
                  ENDFB(G);
              }
            }
          } else {
            PRINTFB(G, FB_Matrix, FB_Warnings)
              "Matrix-Warning: Kabsch eigenvalue(s) negative: falling back on iteration.\n"
              ENDFB(G);
          }
        } else {
          PRINTFB(G, FB_Matrix, FB_Warnings)
            "Matrix-Warning: Kabsch decomposition failed: falling back on iteration.\n"
            ENDFB(G);
        }
      }
    }

    if(!got_kabsch) {

      /* use PyMOL's original iteration algorithm if Kabsch is
         disabled or fails (quite common). */

      /* Primary iteration scheme to determine rotation matrix for molecule 2 */
      double sg, bb, cc;
      int iters, iy, iz, unchanged = 0;
      double sig, gam;
      double tmp;
      double save[3][3];
      int perturbed = false;

      for(a = 0; a < 3; a++) {
        for(b = 0; b < 3; b++) {
          m[a][b] = 0.0;
          save[a][b] = aa[a][b];
        }
        m[a][a] = 1.0;
      }

      iters = 0;
      while(1) {

        /* IX, IY, and IZ rotate 1-2-3, 2-3-1, 3-1-2, etc. */
        iz = (iters + 1) % 3;
        iy = (iz + 1) % 3;
        sig = aa[iz][iy] - aa[iy][iz];
        gam = aa[iy][iy] + aa[iz][iz];

        if(iters >= maxiter) {
          PRINTFB(G, FB_Matrix, FB_Details)
            " Matrix: Warning: no convergence (%1.8f<%1.8f after %d iterations).\n",
            (float) tol, (float) gam, iters ENDFB(G);
          break;
        }

        /* Determine size of off-diagonal element.  If off-diagonals exceed the
           diagonal elements tolerance, perform Jacobi rotation. */
        tmp = sig * sig + gam * gam;
        sg = sqrt1d(tmp);
        if((sg != 0.0F) && (fabs(sig) > (tol * fabs(gam)))) {
          unchanged = 0;
          sg = 1.0F / sg;
          for(a = 0; a < 3; a++) {
            bb = gam * aa[iy][a] + sig * aa[iz][a];
            cc = gam * aa[iz][a] - sig * aa[iy][a];
            aa[iy][a] = bb * sg;
            aa[iz][a] = cc * sg;

            bb = gam * m[iy][a] + sig * m[iz][a];
            cc = gam * m[iz][a] - sig * m[iy][a];
            m[iy][a] = bb * sg;
            m[iz][a] = cc * sg;
          }
        } else {
          unchanged++;
          if(unchanged == 3) {
            double residual = 0.0;
            for(a = 0; a < 3; a++) {
              for(b = 0; b < 3; b++) {
                residual += fabs(aa[a][b] - save[a][b]);
              }
            }
            if(residual > R_SMALL4) {
              /* matrix has changed significantly, so we assume that
                 we found the minimum */
              break;
            } else if(perturbed) {
              /* we ended up back where we started even after perturbing */
              break;
            } else {            /* hmm...no change from start... so displace 90
                                   degrees just to make sure we didn't start out
                                   trapped in precisely the opposite direction */
              for(a = 0; a < 3; a++) {
                bb = aa[iz][a];
                cc = -aa[iy][a];
                aa[iy][a] = bb;
                aa[iz][a] = cc;

                bb = m[iz][a];
                cc = -m[iy][a];
                m[iy][a] = bb;
                m[iz][a] = cc;
              }
              perturbed = true;
              unchanged = 0;
            }
            /* only give up if we've iterated through all three planes */
          }
        }
        iters++;
      }
      recondition33d((double *) (void *) m);
    }
  }

  /* At this point, we should have a converged rotation matrix (M).  Calculate
     the weighted RMS error. */
  {
    const float *vv1 = v1, *vv2 = v2;
    double etmp, tmp;
    double err = 0.0;
    for(c = 0; c < n; c++) {
      etmp = 0.0;
      for(a = 0; a < 3; a++) {
        tmp = m[a][0] * (vv2[0] - t2[0])
          + m[a][1] * (vv2[1] - t2[1])
          + m[a][2] * (vv2[2] - t2[2]);
        tmp = (vv1[a] - t1[a]) - tmp;
        etmp += tmp * tmp;
      }
      if(wt)
        err += wt[c] * etmp;
      else
        err += etmp;
      vv1 += 3;
      vv2 += 3;
    }

    err = err / sumwt;
    err = sqrt1d(err);

    /* NOTE: TTT's are now row-major (to be more like homogenous matrices) */

    if(ttt) {
      ttt[0] = (float) m[0][0];
      ttt[1] = (float) m[1][0];
      ttt[2] = (float) m[2][0];
      ttt[3] = (float) t2[0];
      ttt[4] = (float) m[0][1];
      ttt[5] = (float) m[1][1];
      ttt[6] = (float) m[2][1];
      ttt[7] = (float) t2[1];
      ttt[8] = (float) m[0][2];
      ttt[9] = (float) m[1][2];
      ttt[10] = (float) m[2][2];
      ttt[11] = (float) t2[2];
      ttt[12] = (float) -t1[0];
      ttt[13] = (float) -t1[1];
      ttt[14] = (float) -t1[2];
    }
    /* for compatibility with normal 4x4 matrices */

    if(fabs(err) < R_SMALL4)
      err = 0.0F;

    return ((float) err);
  }
}


/*========================================================================*/

/**
 * @param G Only used for feedback, can be NULL
 * @param[in] a 3x3 matrix
 * @param[out] wr 3x1 eigenvalues (real part of complex number)
 * @param[out] wi 3x1 eigenvalues (imag part of complex number)
 * @param[out] v 3x3 eigenvectors
 * @return Error code
 */
int MatrixEigensolveC33d(PyMOLGlobals * G, const double *a, double *wr, double *wi, double *v)
{
  TNT::Array2D<double> A(3, 3);
  TNT::Array2D<double> V(3, 3);
  TNT::Array1D<double> WR(3);
  TNT::Array1D<double> WI(3);

  // input
  transpose33d33d(a, A[0]);

  JAMA::Eigenvalue<double> E(A);
  E.getRealEigenvalues(WR);
  E.getImagEigenvalues(WI);
  E.getV(V);

  // output
  transpose33d33d(V[0], v);
  copy3d(static_cast<double const*>(WR), wr);
  copy3d(static_cast<double const*>(WI), wi);

  /* NOTE: the returned eigenvectors are stored one per row which is
     is actually the inverse of the normal eigenvalue matrix -
     ----
     IS that because we're actually solving the transpose?
   */

  if (G && Feedback(G, FB_Matrix, FB_Blather)) {
    printf(" Eigensolve: eigenvectors %8.3f %8.3f %8.3f\n", v[0], v[1], v[2]);
    printf(" Eigensolve:              %8.3f %8.3f %8.3f\n", v[3], v[4], v[5]);
    printf(" Eigensolve:              %8.3f %8.3f %8.3f\n", v[6], v[7], v[8]);

    printf(" Eigensolve: eigenvalues  %8.3f %8.3f %8.3f\n", wr[0], wr[1], wr[2]);
    printf(" Eigensolve:              %8.3f %8.3f %8.3f\n", wi[0], wi[1], wi[2]);
  }

  return 0;
}
