/*
 * totvar.c --
 *
 * Implementation of relaxed Total Variation (TV) in 2D or 3D.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2016, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>,
 *                          Loïc Denis <loic.denis@univ-st-etienne.fr>,
 *                          Ferréol Soulez <ferreol.soulez@univ-lyon1.fr>,
 *                          Fabien Momey <fabien.momey@gmail.com>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>

#include "totvar.h"

/* Macro to compute the square of its argument. */
#define SQ(a) ((a)*(a))

/* Macro to check method in options. */
#define METHOD(option, value) (((option)&(7)) == (value))

/* Offsets in multidimensional arrays (dimensions are stored in variables n1,
   n2, etc.). */
#define  OFF1(a1)            (a1)
#define  OFF2(a1,a2)         ((a2)*n1 + (a1))
#define  OFF3(a1,a2,a3)      (((a3)*n2 + (a2))*n1 + (a1))
#define  OFF4(a1,a2,a3,a4)   ((((a4)*n3 + (a3))*n2 + (a2))*n1 + (a1))

/* Derivative of the absolute value of Q*Y w.r.t. Y. */
#define DRV_ABS(y,q) ((y) > ZERO ? (q) : ((y) < ZERO ? -(q) : ZERO))

#define COMPUTE_GRADIENT (RGL_INTEGRATE_GRADIENT | RGL_STORE_GRADIENT)

double rgl_tv2d(const double x[], const long n1, const long n2,
                const double w1, const double w2, const double eps,
                double gx[], unsigned int flags)
{
  const double ZERO = 0.0;
  double x1, x2, x3, x4, y1, y2, y3;
  double p, r, s, fx;
  long i1, i2, j1, j2, j3, j4;
  int compute_gradient;

  compute_gradient = ((flags & COMPUTE_GRADIENT) != 0);
  if (x == NULL || n1 < 0 || n2 < 0 || w1 < 0.0 || w2 < 0.0 || eps < 0.0 ||
      (gx == NULL && compute_gradient)) {
    return -1.0;
  }
  if ((flags & RGL_STORE_GRADIENT) != 0) {
    memset(gx, 0, n1*n2*sizeof(gx[0]));
  }
  fx = 0.0;
  s = eps*eps;

#define  X(a1,a2)   x[OFF2(a1,a2)]
#define GX(a1,a2)  gx[OFF2(a1,a2)]

  if (METHOD(flags, RGL_TOTVAR_ALTERNATIVE)) {

    /*
     * Use "alternative" isotropic definition of the squared norm of the
     * spatial gradient which is the average of the squared norm of the
     * spatial gradient of the bilinear interpolation function.
     *
     * For each position, the cost writes:
     *
     *   F = sqrt(R + S) - EPS
     *
     * with S = EPS*EPS and R is the quadratic norm of the local gradient:
     *
     *   R = Q1*Y1*Y1 + Q2*Y2*Y2 + Q3*Y3*Y3
     *
     * with:
     *
     *   Q1  =  W1/4
     *   Q2  =  W2/4
     *   Q3  =  (W1 + W2)/12
     *
     *   Y1  =  X1 - X2 + X3 - X4  =  (X1 - X4) - (X2 - X3)
     *   Y2  =  X1 + X2 - X3 - X4  =  (X1 - X4) + (X2 - X3)
     *   Y3  =  X1 - X2 - X3 + X4  =  (X1 - X3) - (X2 - X4)
     *
     *             i1
     *        +------->
     *        | X1 X2
     *     i2 | X3 X4
     *        V
     *
     * If all the weights are the same, the gradient becomes:
     *
     *   R  =  Q1*Y1*Y1 + Q2*Y2*Y2 + Q3*Y3*Y3
     *      =  Q1*(Y1*Y1 + Y2*Y2) + Q3*Y3*Y3
     *
     * with:
     *
     *   Q1  =  Q2  =  W/2
     *   Q3  =  W/6
     *   Y1  =  X1 - X4;
     *   Y2  =  X2 - X3;
     *   Y3  =  X1 - X2 - X3 + X4  =  (X1 - X3) - (X2 - X4)
     */

    if (compute_gradient) {
      if (w1 == w2) /* same weights along all directions */ {
        double q1 = w1/2.0;
        double q3 = w1/6.0;
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          x4 = X(0, i2);
          for (i1 = 1; i1 < n1; ++i1) {
            x3 = x4;
            x1 = x2;
            x2 = X(i1, i2-1);
            x4 = X(i1, i2);
            y1 = x1 - x4;
            y2 = x2 - x3;
            y3 = (x1 - x3) - (x2 - x4);
            /* FIXME: it is possible to save 2 multiplications */
            r = q1*(y1*y1 + y2*y2) + q3*y3*y3;
            p = sqrt(r + s);
            fx += p;
            p = 1.0/p;
            y1 *= p*q1;
            y2 *= p*q1;
            y3 *= p*q3;
            GX(i1 - 1, i2 - 1) += y1 + y3; /* deriv. wrt X1 */
            GX(i1    , i2 - 1) += y2 - y3; /* deriv. wrt X2 */
            GX(i1 - 1, i2    ) -= y2 + y3; /* deriv. wrt X3 */
            GX(i1    , i2    ) -= y1 - y3; /* deriv. wrt X4 */
          }
        }
      } else /* not same weights along all directions */ {
        double q1 = w1/4.0;
        double q2 = w2/4.0;
        double q3 = (w1 + w2)/12.0;
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          x4 = X(0, i2);
          for (i1 = 1; i1 < n1; ++i1) {
            x3 = x4;
            x1 = x2;
            x2 = X(i1, i2-1);
            x4 = X(i1, i2);
            y1 = (x1 - x4) - (x2 - x3);
            y2 = (x1 - x4) + (x2 - x3);
            y3 = (x1 - x3) - (x2 - x4);
            r = q1*y1*y1 + q2*y2*y2 + q3*y3*y3;
            p = sqrt(r + s);
            fx += p;
            /* FIXME: it is possible to save 3 multiplications */
            p = 1.0/p;
            y1 *= p*q1;
            y2 *= p*q2;
            y3 *= p*q3;
            GX(i1 - 1, i2 - 1) += y1 + y2 + y3; /* deriv. wrt X1 */
            GX(i1    , i2 - 1) -= y1 - y2 + y3; /* deriv. wrt X2 */
            GX(i1 - 1, i2    ) += y1 - y2 - y3; /* deriv. wrt X3 */
            GX(i1    , i2    ) -= y1 + y2 - y3; /* deriv. wrt X4 */
          }
        }
      }
    } else /* do not compute gradients */ {
      if (w1 == w2) /* same weights along all directions */ {
        double q1 = w1/2.0;
        double q3 = w1/6.0;
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          x4 = X(0, i2);
          for (i1 = 1; i1 < n1; ++i1) {
            x3 = x4;
            x1 = x2;
            x2 = X(i1, i2-1);
            x4 = X(i1, i2);
            y1 = x1 - x4;
            y2 = x2 - x3;
            y3 = (x1 - x3) - (x2 - x4);
            r = q1*(y1*y1 + y2*y2) + q3*y3*y3;
            fx += sqrt(r + s);
          }
        }
      } else /* not same weights along all directions */ {
        double q1 = w1/4.0;
        double q2 = w2/4.0;
        double q3 = (w1 + w2)/12.0;
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          x4 = X(0, i2);
          for (i1 = 1; i1 < n1; ++i1) {
            x3 = x4;
            x1 = x2;
            x2 = X(i1, i2-1);
            x4 = X(i1, i2);
            y1 = (x1 - x4) - (x2 - x3);
            y2 = (x1 - x4) + (x2 - x3);
            y3 = (x1 - x3) - (x2 - x4);
            r = q1*y1*y1 + q2*y2*y2 + q3*y3*y3;
            fx += sqrt(r + s);
          }
        }
      }
    }

  } else if (METHOD(flags, RGL_TOTVAR_FORWARD)) {

    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * forward differences along each direction (standard discretization).
     */

    double y12, y13;

    if (compute_gradient) {
      if (w1 == w2) /* same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1 = x2;
            x2 = X(i1, i2-1);
            x3 = X(i1-1, i2);
            y12 = x1 - x2;
            y13 = x1 - x3;
            r = (SQ(y12) + SQ(y13))*w1;
            p = sqrt(r + s);
            fx += p;
            p = w1/p;
            y12 *= p;
            y13 *= p;
            GX(i1 - 1, i2 - 1) += y12 + y13; /* deriv. wrt X1 */
            GX(i1    , i2 - 1) -= y12;       /* deriv. wrt X2 */
            GX(i1 - 1, i2    ) -= y13;       /* deriv. wrt X3 */
          }
        }
      } else /* not same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1 = x2;
            x2 = X(i1, i2-1);
            x3 = X(i1-1, i2);
            y12 = x1 - x2;
            y13 = x1 - x3;
            r = SQ(y12)*w1 + SQ(y13)*w2;
            p = sqrt(r + s);
            fx += p;
            p = 1.0/p;
            y12 *= p*w1;
            y13 *= p*w2;
            GX(i1 - 1, i2 - 1) += y12 + y13; /* deriv. wrt X1 */
            GX(i1    , i2 - 1) -= y12;       /* deriv. wrt X2 */
            GX(i1 - 1, i2    ) -= y13;       /* deriv. wrt X3 */
          }
        }
      }
    } else /* do not compute gradients */ {
      if (w1 == w2) /* same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1 = x2;
            x2 = X(i1, i2-1);
            x3 = X(i1-1, i2);
            y12 = x1 - x2;
            y13 = x1 - x3;
            r = (SQ(y12) + SQ(y13))*w1;
            fx += sqrt(r + s);
          }
        }
      } else /* not same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2 = X(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1 = x2;
            x2 = X(i1, i2-1);
            x3 = X(i1-1, i2);
            y12 = x1 - x2;
            y13 = x1 - x3;
            r = SQ(y12)*w1 + SQ(y13)*w2;
            fx += sqrt(r + s);
          }
        }
      }
    }

  } else if (METHOD(flags, RGL_TOTVAR_SEPARABLE)) {

    /*
     * The cost is the sum of the (relaxed) absolute values of the finite
     * differences along all dimensions.
     *
     * I do not tried to optimize this part (considering the different cases)
     * since this regularization does not give very good results.  Hence using
     * it is only relevant for comparison.
     */
    double p1, p2, q1, q2, y1, y2, z1, z2;

#define J1 OFF2(i1-1,i2-1)
#define J2 OFF2(i1,  i2-1)
#define J3 OFF2(i1-1,i2)
#define J4 OFF2(i1,  i2)
    if (s != 0.0) {
      /* first row */
      x2 = x[0];
      for (i1 = 1; i1 < n1; ++i1) {
        x1 = x2;
        x2 = x[i1];
        y1 = x1 - x2;
        z1 = w1*y1;
        p1 = sqrt(z1*y1 + s);
        fx += p1;
        if (compute_gradient) {
          z1 /= p1;
          gx[i1 - 1] += z1;
          gx[i1]     -= z1;
        }
      }
      /* other rows */
      for (i2 = 1; i2 < n2; ++i2) {
        i1 = 0;
        x2 = x[J2];
        for (i1 = 1; i1 < n1; ++i1) {
          x1 = x2;
          x2 = x[J2];
          x3 = x[J3];
          y1 = x1 - x2;
          z1 = w1*y1;
          p1 = sqrt(z1*y1 + s);
          y2 = x1 - x3;
          z2 = w2*y2;
          p2 = sqrt(z2*y2 + s);
          fx += p1 + p2;
          if (compute_gradient) {
            z1 /= p1;
            z2 /= p2;
            gx[J1] += z1 + z2;
            gx[J2] -= z1;
            gx[J3] -= z2;
          }
        }
      }
    } else /* s = 0 ==> use absolute value */ {
      /* first row */
      q1 = sqrt(w1);
      q2 = sqrt(w2);
      x2 = x[0];
      for (i1 = 1; i1 < n1; ++i1) {
        x1 = x2;
        x2 = x[i1];
        y1 = x1 - x2;
        p1 = fabs(q1*y1);
        fx += p1;
        if (compute_gradient) {
          y1 = DRV_ABS(y1 , q1);
          gx[i1 - 1] += y1;
          gx[i1]     -= y1;
        }
      }
      /* other rows */
      for (i2 = 1; i2 < n2; ++i2) {
        i1 = 0;
        x2 = x[J2];
       for (i1 = 1; i1 < n1; ++i1) {
          x1 = x2;
          x2 = x[J2];
          x3 = x[J3];
          y1 = x1 - x2;
          p1 = fabs(q1*y1);
          y2 = x1 - x3;
          p2 = fabs(q2*y2);
          fx += p1 + p2;
          if (compute_gradient) {
            y1 = DRV_ABS(y1 , q1);
            y2 = DRV_ABS(y2 , q2);
            gx[J1] += y1 + y2;
            gx[J2] -= y1;
            gx[J3] -= y2;
          }
        }
      }
    }
#undef J1
#undef J2
#undef J3
#undef J4

  } else {

    /*
     * Assume RGL_TOTVAR_ISOTROPIC.
     *
     * The squared norm of the spatial gradient is the sum of the
     * squared differences along all the edges (divided by 2).
     */

    double y12, y34, y13, y24;

    if (w1 == w2) /* same weights along all directions */ {
      double q = w1/2.0;
      if (compute_gradient) {
        double p;
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2 = x[j2];
          x4 = x[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            j1 = j2++;
            j3 = j4++;
            x1 = x2;
            x2 = x[j2];
            x3 = x4;
            x4 = x[j4];
            y12 = x1 - x2;
            y34 = x3 - x4;
            y13 = x1 - x3;
            y24 = x2 - x4;
            r = (SQ(y12) + SQ(y34) + SQ(y13) + SQ(y24))*q;
            p = sqrt(r + s);
            fx += p;
            p = q/p;
            gx[j1] += (y12 + y13)*p;
            gx[j2] -= (y12 - y24)*p;
            gx[j3] += (y34 - y13)*p;
            gx[j4] -= (y34 + y24)*p;
          }
        }
      } else /* do not compute gradient */ {
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2 = x[j2];
          x4 = x[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            x1 = x2;
            x2 = x[++j2];
            x3 = x4;
            x4 = x[++j4];
            y12 = x1 - x2;
            y34 = x3 - x4;
            y13 = x1 - x3;
            y24 = x2 - x4;
            r = (SQ(y12) + SQ(y34) + SQ(y13) + SQ(y24))*q;
            fx += sqrt(r + s);
          }
        }
      }
    } else /* not same weights along all directions */ {
      double q1 = w1/2.0, q2 = w2/2.0;
      if (compute_gradient) {
        double p1, p2;
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2 = x[j2];
          x4 = x[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            j1 = j2++;
            j3 = j4++;
            x1 = x2;
            x2 = x[j2];
            x3 = x4;
            x4 = x[j4];
            y12 = x1 - x2;
            y34 = x3 - x4;
            y13 = x1 - x3;
            y24 = x2 - x4;
            r = (SQ(y12) + SQ(y34))*q1 + (SQ(y13) + SQ(y24))*q2;
            r = sqrt(r + s);
            fx += r;
            r = 1.0/r;
            p1 = r*q1;
            p2 = r*q2;
            y12 *= p1;
            y34 *= p1;
            y13 *= p2;
            y24 *= p2;
            gx[j1] += y12 + y13;
            gx[j2] -= y12 - y24;
            gx[j3] += y34 - y13;
            gx[j4] -= y34 + y24;
          }
        }
      } else /* do not compute gradient */ {
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2 = x[j2];
          x4 = x[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            x1 = x2;
            x2 = x[++j2];
            x3 = x4;
            x4 = x[++j4];
            y12 = x1 - x2;
            y34 = x3 - x4;
            y13 = x1 - x3;
            y24 = x2 - x4;
            r = (SQ(y12) + SQ(y34))*q1 + (SQ(y13) + SQ(y24))*q2;
            fx += sqrt(r + s);
          }
        }
      }
    }
  }

#undef X
#undef GX

  /* Remove the "bias" and make sure the result is non-negative (it can only
     be negative due to rounding errors). */
  fx -= (n1 - 1)*(n2 - 1)*eps;
  if (fx < 0.0) {
    fx = 0.0;
  }
  return fx;
}

double rgl_tv2d_complex(const double xr[], const double xi[],
                       const long n1, const long n2,
                       const double w1, const double w2,
                       const double eps, double gxr[], double gxi[],
                       unsigned int flags)
{
  const double ZERO = 0.0;
  double x1r, x2r, x3r, x4r, y1r, y2r, y3r;
  double x1i, x2i, x3i, x4i, y1i, y2i, y3i;
  double p, r, s, fx;
  long i1, i2, j1, j2, j3, j4;
  int compute_gradient;

  compute_gradient = ((flags & COMPUTE_GRADIENT) != 0);
  if (xr == NULL || xi ==NULL || n1 < 0 || n2 < 0 || w1 < 0.0 || w2 < 0.0 || eps < 0.0 ||
      (gxr == NULL && compute_gradient) || (gxi == NULL && compute_gradient)) {
    return -1.0;
  }
  if ((flags & RGL_STORE_GRADIENT) != 0) {
    memset(gxr, 0, n1*n2*sizeof(gxr[0]));
    memset(gxi, 0, n1*n2*sizeof(gxi[0]));
  }
  fx = 0.0;
  s = eps*eps;

#define  XR(a1,a2)   xr[OFF2(a1,a2)]
#define GXR(a1,a2)  gxr[OFF2(a1,a2)]
#define  XI(a1,a2)   xi[OFF2(a1,a2)]
#define GXI(a1,a2)  gxi[OFF2(a1,a2)]

  if (METHOD(flags, RGL_TOTVAR_FORWARD)) {
    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * forward differences along each direction (standard discretization).
     */

    double y12r, y13r, y12i, y13i;
    
    if (compute_gradient) {
      if (w1 == w2) /* same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2r = XR(0, i2-1);
          x2i = XI(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1r = x2r;           x1i = x2i;
            x2r = XR(i1, i2-1);  x2i = XI(i1, i2-1);
            x3r = XR(i1-1, i2);  x3i = XI(i1-1, i2);
            y12r = x1r - x2r;    y12i = x1i - x2i;
            y13r = x1r - x3r;    y13i = x1i - x3i;
            r = (SQ(y12r) + SQ(y13r) + SQ(y12i) + SQ(y13i))*w1;
            p = sqrt(r + s);
            fx += p;
            p = w1/p;
            y12r *= p;
            y13r *= p;
            y12i *= p;
            y13i *= p;
            GXR(i1 - 1, i2 - 1) += y12r + y13r; /* deriv. wrt X1 */
            GXR(i1    , i2 - 1) -= y12r;       /* deriv. wrt X2 */
            GXR(i1 - 1, i2    ) -= y13r;       /* deriv. wrt X3 */
            GXI(i1 - 1, i2 - 1) += y12i + y13i; /* deriv. wrt X1 */
            GXI(i1    , i2 - 1) -= y12i;       /* deriv. wrt X2 */
            GXI(i1 - 1, i2    ) -= y13i;       /* deriv. wrt X3 */
          }
        }
      } else /* not same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2r = XR(0, i2-1);
          x2i = XI(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1r = x2r;           x1i = x2i;
            x2r = XR(i1, i2-1);  x2i = XI(i1, i2-1);
            x3r = XR(i1-1, i2);  x3i = XI(i1-1, i2);
            y12r = x1r - x2r;    y12i = x1i - x2i;
            y13r = x1r - x3r;    y13i = x1i - x3i;
            r = SQ(y12r)*w1 + SQ(y13r)*w2 + SQ(y12i)*w1 + SQ(y13i)*w2;
            p = sqrt(r + s);
            fx += p;
            p = 1.0/p;
            y12r *= p*w1;
            y13r *= p*w2;
            y12i *= p*w1;
            y13i *= p*w2;
            GXR(i1 - 1, i2 - 1) += y12r + y13r; /* deriv. wrt X1 */
            GXR(i1    , i2 - 1) -= y12r;       /* deriv. wrt X2 */
            GXR(i1 - 1, i2    ) -= y13r;       /* deriv. wrt X3 */
            GXI(i1 - 1, i2 - 1) += y12i + y13i; /* deriv. wrt X1 */
            GXI(i1    , i2 - 1) -= y12i;       /* deriv. wrt X2 */
            GXI(i1 - 1, i2    ) -= y13i;       /* deriv. wrt X3 */
          }
        }
      }
    } else /* do not compute gradients */ {
      if (w1 == w2) /* same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2r = XR(0, i2-1);
          x2i = XI(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1r = x2r;           x1i = x2i;
            x2r = XR(i1, i2-1);  x2i = XI(i1, i2-1);
            x3r = XR(i1-1, i2);  x3i = XI(i1-1, i2);
            y12r = x1r - x2r;    y12i = x1i - x2i;
            y13r = x1r - x3r;    y13i = x1i - x3i;
            r = (SQ(y12r) + SQ(y13r) + SQ(y12i) + SQ(y13i))*w1;
            fx += sqrt(r + s);
          }
        }
      } else /* not same weights along all directions */ {
        for (i2 = 1; i2 < n2; ++i2) {
          x2r = XR(0, i2-1);
          x2i = XI(0, i2-1);
          for (i1 = 1; i1 < n1; ++i1) {
            x1r = x2r;           x1i = x2i;
            x2r = XR(i1, i2-1);  x2i = XI(i1, i2-1);
            x3r = XR(i1-1, i2);  x3i = XI(i1-1, i2);
            y12r = x1r - x2r;    y12i = x1i - x2i;
            y13r = x1r - x3r;    y13i = x1i - x3i;
            r = SQ(y12r)*w1 + SQ(y13r)*w2 + SQ(y12i)*w1 + SQ(y13i)*w2;
            fx += sqrt(r + s);
          }
        }
      }
    }
  
  } else if (METHOD(flags, RGL_TOTVAR_ISOTROPIC)) {
    /*
     * Assume RGL_TOTVAR_ISOTROPIC.
     *
     * The squared norm of the spatial gradient is the sum of the
     * squared differences along all the edges (divided by 2).
     */

    double y12r, y34r, y13r, y24r;
    double y12i, y34i, y13i, y24i;

    if (w1 == w2) /* same weights along all directions */ {
      double q = w1/2.0;
      if (compute_gradient) {
        double p;
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2r = xr[j2]; x2i = xi[j2]; 
          x4r = xr[j4]; x4i = xi[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            j1 = j2++;
            j3 = j4++;
            x1r = x2r;     x1i = x2i;
            x2r = xr[j2];  x2i = xi[j2];
            x3r = x4r;     x3i = x4i;
            x4r = xr[j4];  x4i = xi[j4];
            y12r = x1r - x2r;  y12i = x1i - x2i;
            y34r = x3r - x4r;  y34i = x3i - x4i;
            y13r = x1r - x3r;  y13i = x1i - x3i;
            y24r = x2r - x4r;  y24i = x2i - x4i;
            r = (SQ(y12r) + SQ(y34r) + SQ(y13r) + SQ(y24r) +
                 SQ(y12i) + SQ(y34i) + SQ(y13i) + SQ(y24i))*q;
            p = sqrt(r + s);
            fx += p;
            p = q/p;
            gxr[j1] += (y12r + y13r)*p;
            gxr[j2] -= (y12r - y24r)*p;
            gxr[j3] += (y34r - y13r)*p;
            gxr[j4] -= (y34r + y24r)*p;
            gxi[j1] += (y12i + y13i)*p;
            gxi[j2] -= (y12i - y24i)*p;
            gxi[j3] += (y34i - y13i)*p;
            gxi[j4] -= (y34i + y24i)*p;
          }
        }
      } else /* do not compute gradient */ {
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2r = xr[j2]; x2i = xi[j2];
          x4r = xr[j4]; x4i = xi[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            j2++;
            j4++;
            x1r = x2r;       x1i = x2i;
            x2r = xr[j2];  x2i = xi[j2];
            x3r = x4r;       x3i = x4i;
            x4r = xr[j4];  x4i = xi[j4];
            y12r = x1r - x2r;  y12i = x1i - x2i;
            y34r = x3r - x4r;  y34i = x3i - x4i;
            y13r = x1r - x3r;  y13i = x1i - x3i;
            y24r = x2r - x4r;  y24i = x2i - x4i;
            r = (SQ(y12r) + SQ(y34r) + SQ(y13r) + SQ(y24r) +
                 SQ(y12i) + SQ(y34i) + SQ(y13i) + SQ(y24i))*q;
            fx += sqrt(r + s);
          }
        }
      }
    } else /* not same weights along all directions */ {
      double q1 = w1/2.0, q2 = w2/2.0;
      if (compute_gradient) {
        double p1, p2;
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2r = xr[j2];  x2i = xi[j2];
          x4r = xr[j4];  x4i = xi[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            j1 = j2++;
            j3 = j4++;
            x1r = x2r;     x1i = x2i;
            x2r = xr[j2];  x2i = xi[j2];
            x3r = x4r;     x3i = x4i;
            x4r = xr[j4];  x4i = xi[j4];
            y12r = x1r - x2r;  y12i = x1i - x2i;
            y34r = x3r - x4r;  y34i = x3i - x4i;
            y13r = x1r - x3r;  y13i = x1i - x3i;
            y24r = x2r - x4r;  y24i = x2i - x4i;
            r = (SQ(y12r) + SQ(y34r))*q1 + (SQ(y13r) + SQ(y24r))*q2 +
              (SQ(y12i) + SQ(y34i))*q1 + (SQ(y13i) + SQ(y24i))*q2;
            r = sqrt(r + s);
            fx += r;
            r = 1.0/r;
            p1 = r*q1;
            p2 = r*q2;
            y12r *= p1; y12i *= p1;
            y34r *= p1; y34i *= p1;
            y13r *= p2; y13i *= p2;
            y24r *= p2; y24i *= p2;
            gxr[j1] += y12r + y13r;
            gxr[j2] -= y12r - y24r;
            gxr[j3] += y34r - y13r;
            gxr[j4] -= y34r + y24r;
            gxi[j1] += y12i + y13i;
            gxi[j2] -= y12i - y24i;
            gxi[j3] += y34i - y13i;
            gxi[j4] -= y34i + y24i;
          }
        }
      } else /* do not compute gradient */ {
        for (i2 = 1; i2 < n2; ++i2) {
          j2 = OFF2(0, i2 - 1);
          j4 = OFF2(0, i2);
          x2r = xr[j2]; x2i = xi[j2];
          x4r = xr[j4]; x4i = xi[j4];
          for (i1 = 1; i1 < n1; ++i1) {
            j2++;
            j4++;
            x1r = x2r;     x1i = x2i;
            x2r = xr[j2];  x2i = xi[j2];
            x3r = x4r;     x3i = x4i;
            x4r = xr[j4];  x4i = xi[j4];
            y12r = x1r - x2r;  y12i = x1i - x2i;
            y34r = x3r - x4r;  y34i = x3i - x4i;
            y13r = x1r - x3r;  y13i = x1i - x3i;
            y24r = x2r - x4r;  y24i = x2i - x4i;
            r = (SQ(y12r) + SQ(y34r))*q1 + (SQ(y13r) + SQ(y24r))*q2 +
              (SQ(y12i) + SQ(y34i))*q1 + (SQ(y13i) + SQ(y24i))*q2;
            fx += sqrt(r + s);
          }
        }
      }
    }
  }

#undef X
#undef GX

  /* Remove the "bias" and make sure the result is non-negative (it can only
     be negative due to rounding errors). */
  fx -= (n1 - 1)*(n2 - 1)*eps;
  if (fx < 0.0) {
    fx = 0.0;
  }
  return fx;
}

/*
 * NOTATIONS FOR 3-D VOLUME X(i1,i2,i3)
 *
 *                 i3  i2
 *                  | /
 *                  |/              X1 = X(i1-1,i2-1,i3-1)
 *        X7--------X8---> i1       X2 = X(i1  ,i2-1,i3-1)
 *       /:        /|               X3 = X(i1-1,i2  ,i3-1)
 *      / :       / |               X4 = X(i1  ,i2  ,i3-1)
 *     X5--------X6 |               X5 = X(i1-1,i2-1,i3  )
 *     |  X3.....|..X4              X6 = X(i1  ,i2-1,i3  )
 *     | '       | /                X7 = X(i1-1,i2  ,i3  )
 *     |'        |/                 X8 = X(i1  ,i2  ,i3  )
 *     X1--------X2
 *
 */
double rgl_tv3d(const double x[],
                const long n1, const long n2, const long n3,
                const double w1, const double w2, const double w3,
                const double eps, double gx[], unsigned int flags)
{
  double x1, x2, x3, x4, x5, x6, x7, x8;
  double r, s, fx;
  double y12, y34, y56, y78;
  double y13, y24, y57, y68;
  double y15, y26, y37, y48;
  long i1, i2, i3;
  long j1, j2, j3, j4, j5, j6, j7, j8;
  int compute_gradient;

  compute_gradient = ((flags & COMPUTE_GRADIENT) != 0);
  if (x == NULL || n1 < 0 || n2 < 0 || n3 < 0 ||
      w1 < 0.0 || w2 < 0.0 || w3 < 0.0 || eps < 0.0 ||
      (gx == NULL && compute_gradient)) {
    return -1.0;
  }
  if ((flags & RGL_STORE_GRADIENT) != 0) {
    memset(gx, 0, n1*n2*n3*sizeof(gx[0]));
  }
  fx = 0.0;
  s = eps*eps;

  if (METHOD(flags, RGL_TOTVAR_FORWARD)) {

    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * forward differences along each direction (standard discretization).
     */

    if (w1 == w2 && w2 == w3) /* same weights along all directions */ {
      double p, q = w1;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2 = x[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at 4 corners of the cube */
              j1 = j2++;
              x1 = x2;
              x2 = x[j2];
              x3 = x[j3];
              x5 = x[j5];
              /* Compute the differences along the 3 directions: */
              y12 = x1 - x2;
              y13 = x1 - x3;
              y15 = x1 - x5;
              /* Compute the cost and integrate its gradient. */
              r = (SQ(y12) + SQ(y13) + SQ(y15))*q;
              p = sqrt(r + s);
              fx += p;
              p = q/p;
              gx[j1] += (y12 + y13 + y15)*p;
              gx[j2] -= y12*p;
              gx[j3] -= y13*p;
              gx[j5] -= y15*p;
              ++j3;
              ++j5;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2 = x[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at 4 corners of the cube */
              x1 = x2;
              x2 = x[++j2];
              x3 = x[j3++];
              x5 = x[j5++];
              /* Compute the differences along the 3 directions: */
              y12 = x1 - x2;
              y13 = x1 - x3;
              y15 = x1 - x5;
              /* Compute the cost and integrate its gradient. */
              r = (SQ(y12) + SQ(y13) + SQ(y15))*q;
              fx += sqrt(r + s);
            }
          }
        }
      }
    } else /* not same weights along all directions */ {
      double q1 = w1;
      double q2 = w2;
      double q3 = w3;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2 = x[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at 4 corners of the cube */
              j1 = j2++;
              x1 = x2;
              x2 = x[j2];
              x3 = x[j3];
              x5 = x[j5];
              /* Compute the differences along the 3 directions: */
              y12 = x1 - x2;
              y13 = x1 - x3;
              y15 = x1 - x5;
              /* Compute the cost and integrate its gradient. */
              r =  SQ(y12)*q1 + SQ(y13)*q2 + SQ(y15)*q3;
              r = sqrt(r + s);
              fx += r;
              r = 1.0/r;
              y12 *= r*q1;
              y13 *= r*q2;
              y15 *= r*q3;
              gx[j1] += y12 + y13 + y15;
              gx[j2] -= y12;
              gx[j3] -= y13;
              gx[j5] -= y15;
              ++j3;
              ++j5;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2 = x[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube */
              x1 = x2;
              x2 = x[++j2];
              x3 = x[j3++];
              x5 = x[j5++];
              /* Compute the differences along the 3 directions: */
              y12 = x1 - x2;
              y13 = x1 - x3;
              y15 = x1 - x5;
              /* Compute the cost and integrate its gradient. */
              r = SQ(y12)*q1 + SQ(y13)*q2 + SQ(y15)*q3;
              fx += sqrt(r + s);
            }
          }
        }
      }
    }

  } else /* same as RGL_TOTVAR_ISOTROPIC */ {

    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * differences along all egdes of the 2x2x2 cube divided by 4 (isotropic
     * definition).
     */

    if (w1 == w2 && w2 == w3) /* same weights along all directions */ {
      double p, q = w1/4.0;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2 = x[j2];
            x4 = x[j4];
            x6 = x[j6];
            x8 = x[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube */
              j1 = j2++; x1 = x2; x2 = x[j2];
              j3 = j4++; x3 = x4; x4 = x[j4];
              j5 = j6++; x5 = x6; x6 = x[j6];
              j7 = j8++; x7 = x8; x8 = x[j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12 = x1 - x2;
              y34 = x3 - x4;
              y56 = x5 - x6;
              y78 = x7 - x8;
              /*  - along 2nd dim: */
              y13 = x1 - x3;
              y24 = x2 - x4;
              y57 = x5 - x7;
              y68 = x6 - x8;
              /*  - along 3rd dim: */
              y15 = x1 - x5;
              y26 = x2 - x6;
              y37 = x3 - x7;
              y48 = x4 - x8;
              /* Compute the cost and integrate its gradient. */
              r = (SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78) +
                   SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68) +
                   SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48))*q;
              p = sqrt(r + s);
              fx += p;
              p = q/p;
              gx[j1] += (y12 + y13 + y15)*p;
              gx[j2] -= (y12 - y24 - y26)*p;
              gx[j3] += (y34 - y13 + y37)*p;
              gx[j4] -= (y34 + y24 - y48)*p;
              gx[j5] += (y56 + y57 - y15)*p;
              gx[j6] -= (y56 - y68 + y26)*p;
              gx[j7] += (y78 - y57 - y37)*p;
              gx[j8] -= (y78 + y68 + y48)*p;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2 = x[j2];
            x4 = x[j4];
            x6 = x[j6];
            x8 = x[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube */
              x1 = x2; x2 = x[++j2];
              x3 = x4; x4 = x[++j4];
              x5 = x6; x6 = x[++j6];
              x7 = x8; x8 = x[++j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12 = x1 - x2;
              y34 = x3 - x4;
              y56 = x5 - x6;
              y78 = x7 - x8;
              /*  - along 2nd dim: */
              y13 = x1 - x3;
              y24 = x2 - x4;
              y57 = x5 - x7;
              y68 = x6 - x8;
              /*  - along 3rd dim: */
              y15 = x1 - x5;
              y26 = x2 - x6;
              y37 = x3 - x7;
              y48 = x4 - x8;
              /* Compute the cost. */
              r = (SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78) +
                   SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68) +
                   SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48))*q;
              fx += sqrt(r + s);
            }
          }
        }
      }
    } else /* not same weights along all directions */ {
      double p1, q1 = w1/4.0;
      double p2, q2 = w2/4.0;
      double p3, q3 = w3/4.0;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2 = x[j2];
            x4 = x[j4];
            x6 = x[j6];
            x8 = x[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube. */
              j1 = j2++; x1 = x2; x2 = x[j2];
              j3 = j4++; x3 = x4; x4 = x[j4];
              j5 = j6++; x5 = x6; x6 = x[j6];
              j7 = j8++; x7 = x8; x8 = x[j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12 = x1 - x2;
              y34 = x3 - x4;
              y56 = x5 - x6;
              y78 = x7 - x8;
              /*  - along 2nd dim: */
              y13 = x1 - x3;
              y24 = x2 - x4;
              y57 = x5 - x7;
              y68 = x6 - x8;
              /*  - along 3rd dim: */
              y15 = x1 - x5;
              y26 = x2 - x6;
              y37 = x3 - x7;
              y48 = x4 - x8;
              /* Compute the cost and integrate its gradient. */
              r = ((SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78))*q1 +
                   (SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68))*q2 +
                   (SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48))*q3);
              r = sqrt(r + s);
              fx += r;
              r = 1.0/r;
              p1 = r*q1;
              p2 = r*q2;
              p3 = r*q3;
              y12 *= p1;
              y34 *= p1;
              y56 *= p1;
              y78 *= p1;
              y13 *= p2;
              y24 *= p2;
              y57 *= p2;
              y68 *= p2;
              y15 *= p3;
              y26 *= p3;
              y37 *= p3;
              y48 *= p3;
              gx[j1] += y12 + y13 + y15;
              gx[j2] -= y12 - y24 - y26;
              gx[j3] += y34 - y13 + y37;
              gx[j4] -= y34 + y24 - y48;
              gx[j5] += y56 + y57 - y15;
              gx[j6] -= y56 - y68 + y26;
              gx[j7] += y78 - y57 - y37;
              gx[j8] -= y78 + y68 + y48;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2 = x[j2];
            x4 = x[j4];
            x6 = x[j6];
            x8 = x[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube. */
              j1 = j2++; x1 = x2; x2 = x[j2];
              j3 = j4++; x3 = x4; x4 = x[j4];
              j5 = j6++; x5 = x6; x6 = x[j6];
              j7 = j8++; x7 = x8; x8 = x[j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12 = x1 - x2;
              y34 = x3 - x4;
              y56 = x5 - x6;
              y78 = x7 - x8;
              /*  - along 2nd dim: */
              y13 = x1 - x3;
              y24 = x2 - x4;
              y57 = x5 - x7;
              y68 = x6 - x8;
              /*  - along 3rd dim: */
              y15 = x1 - x5;
              y26 = x2 - x6;
              y37 = x3 - x7;
              y48 = x4 - x8;
              /* Compute the cost and integrate its gradient. */
              r = ((SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78))*q1 +
                   (SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68))*q2 +
                   (SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48))*q3);
              fx += sqrt(r + s);
            }
          }
        }
      }
    }

  }

  /* Remove the "bias" and make sure the result is non-negative (it can only
     be negative due to rounding errors). */
  fx -= (n1 - 1)*(n2 - 1)*(n3 - 1)*eps;
  if (fx < 0.0) {
    fx = 0.0;
  }
  return fx;
}

/*
 * NOTATIONS FOR 3-D VOLUME X(i1,i2,i3)
 *
 *                 i3  i2
 *                  | /
 *                  |/              X1 = X(i1-1,i2-1,i3-1)
 *        X7--------X8---> i1       X2 = X(i1  ,i2-1,i3-1)
 *       /:        /|               X3 = X(i1-1,i2  ,i3-1)
 *      / :       / |               X4 = X(i1  ,i2  ,i3-1)
 *     X5--------X6 |               X5 = X(i1-1,i2-1,i3  )
 *     |  X3.....|..X4              X6 = X(i1  ,i2-1,i3  )
 *     | '       | /                X7 = X(i1-1,i2  ,i3  )
 *     |'        |/                 X8 = X(i1  ,i2  ,i3  )
 *     X1--------X2
 *
 */
double rgl_tv3d_complex(const double xr[], const double xi[],
                const long n1, const long n2, const long n3,
                const double w1, const double w2, const double w3,
                const double eps, double gxr[], double gxi[], unsigned int flags)
{
  double x1r, x2r, x3r, x4r, x5r, x6r, x7r, x8r;
  double x1i, x2i, x3i, x4i, x5i, x6i, x7i, x8i;
  double r, s, fx;
  double y12r, y34r, y56r, y78r;
  double y13r, y24r, y57r, y68r;
  double y15r, y26r, y37r, y48r;
  double y12i, y34i, y56i, y78i;
  double y13i, y24i, y57i, y68i;
  double y15i, y26i, y37i, y48i;
  long i1, i2, i3;
  long j1, j2, j3, j4, j5, j6, j7, j8;
  int compute_gradient;

  compute_gradient = ((flags & COMPUTE_GRADIENT) != 0);
  if (xr == NULL || xi == NULL || n1 < 0 || n2 < 0 || n3 < 0 ||
      w1 < 0.0 || w2 < 0.0 || w3 < 0.0 || eps < 0.0 ||
      (gxr == NULL && compute_gradient) || (gxi == NULL && compute_gradient)) {
    return -1.0;
  }
  if ((flags & RGL_STORE_GRADIENT) != 0) {
    memset(gxr, 0, n1*n2*n3*sizeof(gxr[0]));
    memset(gxi, 0, n1*n2*n3*sizeof(gxi[0]));
  }
  fx = 0.0;
  s = eps*eps;

  if (METHOD(flags, RGL_TOTVAR_FORWARD)) {

    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * forward differences along each direction (standard discretization).
     */

    if (w1 == w2 && w2 == w3) /* same weights along all directions */ {
      double p, q = w1;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2r = xr[j2];
            x2i = xi[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at 4 corners of the cube */
              j1 = j2++;
              x1r = x2r;  x1i = x2i;
              x2r = xr[j2];  x2i = xi[j2];
              x3r = xr[j3];  x3i = xi[j3];
              x5r = xr[j5];  x5i = xi[j5];
              /* Compute the differences along the 3 directions: */
              y12r = x1r - x2r;  y12i = x1i - x2i;
              y13r = x1r - x3r;  y13i = x1i - x3i;
              y15r = x1r - x5r;  y15i = x1i - x5i;
              /* Compute the cost and integrate its gradient. */
              r = (SQ(y12r) + SQ(y13r) + SQ(y15r) + SQ(y12i) + SQ(y13i) + SQ(y15i))*q;
              p = sqrt(r + s);
              fx += p;
              p = q/p;
              gxr[j1] += (y12r + y13r + y15r)*p;
              gxr[j2] -= y12r*p;
              gxr[j3] -= y13r*p;
              gxr[j5] -= y15r*p;
              gxi[j1] += (y12i + y13i + y15i)*p;
              gxi[j2] -= y12i*p;
              gxi[j3] -= y13i*p;
              gxi[j5] -= y15i*p;
              ++j3;
              ++j5;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2r = xr[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at 4 corners of the cube */
              x1r = x2r;  x1i = x2i;
              x2r = xr[++j2];  x2i = xi[++j2];
              x3r = xr[j3++];  x3i = xi[j3++];
              x5r = xr[j5++];  x5i = xi[j5++];
              /* Compute the differences along the 3 directions: */
              y12r = x1r - x2r;  y12i = x1i - x2i;
              y13r = x1r - x3r;  y13i = x1i - x3i;
              y15r = x1r - x5r;  y15i = x1i - x5i;
              /* Compute the cost and integrate its gradient. */
              r = (SQ(y12r) + SQ(y13r) + SQ(y15r) + SQ(y12i) + SQ(y13i) + SQ(y15i))*q;
              fx += sqrt(r + s);
            }
          }
        }
      }
    } else /* not same weights along all directions */ {
      double q1 = w1;
      double q2 = w2;
      double q3 = w3;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2r = xr[j2];
            x2i = xi[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at 4 corners of the cube */
              j1 = j2++;
              x1r = x2r;  x1i = x2i;
              x2r = xr[j2];  x2i = xi[j2];
              x3r = xr[j3];  x3i = xi[j3];
              x5r = xr[j5];  x5i = xi[j5];
              /* Compute the differences along the 3 directions: */
              y12r = x1r - x2r;  y12i = x1i - x2i;
              y13r = x1r - x3r;  y13i = x1i - x3i;
              y15r = x1r - x5r;  y15i = x1i - x5i;
              /* Compute the cost and integrate its gradient. */
              r =  SQ(y12r)*q1 + SQ(y13r)*q2 + SQ(y15r)*q3 + SQ(y12i)*q1 + SQ(y13i)*q2 + SQ(y15i)*q3;
              r = sqrt(r + s);
              fx += r;
              r = 1.0/r;
              y12r *= r*q1;  y12i *= r*q1;
              y13r *= r*q2;  y13i *= r*q2;
              y15r *= r*q3;  y15i *= r*q3;
              gxr[j1] += y12r + y13r + y15r;
              gxr[j2] -= y12r;
              gxr[j3] -= y13r;
              gxr[j5] -= y15r;
              gxi[j1] += y12i + y13i + y15i;
              gxi[j2] -= y12i;
              gxi[j3] -= y13i;
              gxi[j5] -= y15i;
              ++j3;
              ++j5;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j2 = OFF3(0, i2-1, i3-1);
            j3 = OFF3(0, i2,   i3-1);
            j5 = OFF3(0, i2-1, i3);
            x2r = xr[j2];
            x2i = xi[j2];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube */
              x1r = x2r;       x1i = x2i;
              x2r = xr[++j2];  x2i = xi[++j2];
              x3r = xr[j3++];  x3i = xi[j3++];
              x5r = xr[j5++];  x5i = xi[j5++];
              /* Compute the differences along the 3 directions: */
              y12r = x1r - x2r;  y12i = x1i - x2i;
              y13r = x1r - x3r;  y13i = x1i - x3i;
              y15r = x1r - x5r;  y15i = x1i - x5i;
              /* Compute the cost and integrate its gradient. */
              r = SQ(y12r)*q1 + SQ(y13r)*q2 + SQ(y15r)*q3 + SQ(y12i)*q1 + SQ(y13i)*q2 + SQ(y15i)*q3;
              fx += sqrt(r + s);
            }
          }
        }
      }
    }

  } else if (METHOD(flags, RGL_TOTVAR_ISOTROPIC)) {

    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * differences along all egdes of the 2x2x2 cube divided by 4 (isotropic
     * definition).
     */

    if (w1 == w2 && w2 == w3) /* same weights along all directions */ {
      double p, q = w1/4.0;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2r = xr[j2]; x2i = xi[j2]; 
            x4r = xr[j4]; x4i = xi[j4];
            x6r = xr[j6]; x6i = xi[j6];
            x8r = xr[j8]; x8i = xi[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube */
              j1 = j2++; x1r = x2r; x2r = xr[j2]; x1i = x2i; x2i = xi[j2];
              j3 = j4++; x3r = x4r; x4r = xr[j4]; x3i = x4i; x4i = xi[j4];
              j5 = j6++; x5r = x6r; x6r = xr[j6]; x5i = x6i; x6i = xi[j6];
              j7 = j8++; x7r = x8r; x8r = xr[j8]; x7i = x8i; x8i = xi[j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12r = x1r - x2r; y12i = x1i - x2i;
              y34r = x3r - x4r; y34i = x3i - x4i;
              y56r = x5r - x6r; y56i = x5i - x6i;
              y78r = x7r - x8r; y78i = x7i - x8i;
              /*  - along 2nd dim: */
              y13r = x1r - x3r; y13i = x1i - x3i;
              y24r = x2r - x4r; y24i = x2i - x4i;
              y57r = x5r - x7r; y57i = x5i - x7i;
              y68r = x6r - x8r; y68i = x6i - x8i;
              /*  - along 3rd dim: */
              y15r = x1r - x5r; y15i = x1i - x5i;
              y26r = x2r - x6r; y26i = x2i - x6i;
              y37r = x3r - x7r; y37i = x3i - x7i;
              y48r = x4r - x8r; y48i = x4i - x8i;
              /* Compute the cost and integrate its gradient. */
              r = (SQ(y12r) + SQ(y34r) + SQ(y56r) + SQ(y78r) +
                   SQ(y13r) + SQ(y24r) + SQ(y57r) + SQ(y68r) +
                   SQ(y15r) + SQ(y26r) + SQ(y37r) + SQ(y48r) +
                   SQ(y12i) + SQ(y34i) + SQ(y56i) + SQ(y78i) +
                   SQ(y13i) + SQ(y24i) + SQ(y57i) + SQ(y68i) +
                   SQ(y15i) + SQ(y26i) + SQ(y37i) + SQ(y48i))*q;
              p = sqrt(r + s);
              fx += p;
              p = q/p;
              gxr[j1] += (y12r + y13r + y15r)*p;
              gxr[j2] -= (y12r - y24r - y26r)*p;
              gxr[j3] += (y34r - y13r + y37r)*p;
              gxr[j4] -= (y34r + y24r - y48r)*p;
              gxr[j5] += (y56r + y57r - y15r)*p;
              gxr[j6] -= (y56r - y68r + y26r)*p;
              gxr[j7] += (y78r - y57r - y37r)*p;
              gxr[j8] -= (y78r + y68r + y48r)*p;
              gxi[j1] += (y12i + y13i + y15i)*p;
              gxi[j2] -= (y12i - y24i - y26i)*p;
              gxi[j3] += (y34i - y13i + y37i)*p;
              gxi[j4] -= (y34i + y24i - y48i)*p;
              gxi[j5] += (y56i + y57i - y15i)*p;
              gxi[j6] -= (y56i - y68i + y26i)*p;
              gxi[j7] += (y78i - y57i - y37i)*p;
              gxi[j8] -= (y78i + y68i + y48i)*p;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2r = xr[j2]; x2i = xi[j2]; 
            x4r = xr[j4]; x4i = xi[j4];
            x6r = xr[j6]; x6i = xi[j6];
            x8r = xr[j8]; x8i = xi[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube */
              x1r = x2r; x2r = xr[++j2]; x1i = x2i; x2i = xi[++j2];
              x3r = x4r; x4r = xr[++j4]; x3i = x4i; x4i = xi[++j4];
              x5r = x6r; x6r = xr[++j6]; x5i = x6i; x6i = xi[++j6];
              x7r = x8r; x8r = xr[++j8]; x7i = x8i; x8i = xi[++j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12r = x1r - x2r; y12i = x1i - x2i;
              y34r = x3r - x4r; y34i = x3i - x4i;
              y56r = x5r - x6r; y56i = x5i - x6i;
              y78r = x7r - x8r; y78i = x7i - x8i;
              /*  - along 2nd dim: */
              y13r = x1r - x3r; y13i = x1i - x3i;
              y24r = x2r - x4r; y24i = x2i - x4i;
              y57r = x5r - x7r; y57i = x5i - x7i;
              y68r = x6r - x8r; y68i = x6i - x8i;
              /*  - along 3rd dim: */
              y15r = x1r - x5r; y15i = x1i - x5i;
              y26r = x2r - x6r; y26i = x2i - x6i;
              y37r = x3r - x7r; y37i = x3i - x7i;
              y48r = x4r - x8r; y48i = x4i - x8i;
              /* Compute the cost. */
              r = (SQ(y12r) + SQ(y34r) + SQ(y56r) + SQ(y78r) +
                   SQ(y13r) + SQ(y24r) + SQ(y57r) + SQ(y68r) +
                   SQ(y15r) + SQ(y26r) + SQ(y37r) + SQ(y48r) +
                   SQ(y12i) + SQ(y34i) + SQ(y56i) + SQ(y78i) +
                   SQ(y13i) + SQ(y24i) + SQ(y57i) + SQ(y68i) +
                   SQ(y15i) + SQ(y26i) + SQ(y37i) + SQ(y48i))*q;
              fx += sqrt(r + s);
            }
          }
        }
      }
    } else /* not same weights along all directions */ {
      double p1, q1 = w1/4.0;
      double p2, q2 = w2/4.0;
      double p3, q3 = w3/4.0;
      if (compute_gradient) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2r = xr[j2]; x2i = xi[j2]; 
            x4r = xr[j4]; x4i = xi[j4];
            x6r = xr[j6]; x6i = xi[j6];
            x8r = xr[j8]; x8i = xi[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube. */
              j1 = j2++; x1r = x2r; x2r = xr[j2]; x1i = x2i; x2i = xi[j2];
              j3 = j4++; x3r = x4r; x4r = xr[j4]; x3i = x4i; x4i = xi[j4];
              j5 = j6++; x5r = x6r; x6r = xr[j6]; x5i = x6i; x6i = xi[j6];
              j7 = j8++; x7r = x8r; x8r = xr[j8]; x7i = x8i; x8i = xi[j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12r = x1r - x2r; y12i = x1i - x2i;
              y34r = x3r - x4r; y34i = x3i - x4i;
              y56r = x5r - x6r; y56i = x5i - x6i;
              y78r = x7r - x8r; y78i = x7i - x8i;
              /*  - along 2nd dim: */
              y13r = x1r - x3r; y13i = x1i - x3i;
              y24r = x2r - x4r; y24i = x2i - x4i;
              y57r = x5r - x7r; y57i = x5i - x7i;
              y68r = x6r - x8r; y68i = x6i - x8i;
              /*  - along 3rd dim: */
              y15r = x1r - x5r; y15i = x1i - x5i;
              y26r = x2r - x6r; y26i = x2i - x6i;
              y37r = x3r - x7r; y37i = x3i - x7i;
              y48r = x4r - x8r; y48i = x4i - x8i;
              /* Compute the cost and integrate its gradient. */
              r = ((SQ(y12r) + SQ(y34r) + SQ(y56r) + SQ(y78r))*q1 +
                   (SQ(y13r) + SQ(y24r) + SQ(y57r) + SQ(y68r))*q2 +
                   (SQ(y15r) + SQ(y26r) + SQ(y37r) + SQ(y48r))*q3 +
                   (SQ(y12i) + SQ(y34i) + SQ(y56i) + SQ(y78i))*q1 +
                   (SQ(y13i) + SQ(y24i) + SQ(y57i) + SQ(y68i))*q2 +
                   (SQ(y15i) + SQ(y26i) + SQ(y37i) + SQ(y48i))*q3);
              r = sqrt(r + s);
              fx += r;
              r = 1.0/r;
              p1 = r*q1;
              p2 = r*q2;
              p3 = r*q3;
              y12r *= p1; y12i *= p1;
              y34r *= p1; y34i *= p1;
              y56r *= p1; y56i *= p1;
              y78r *= p1; y78i *= p1;
              y13r *= p2; y13i *= p2;
              y24r *= p2; y24i *= p2;
              y57r *= p2; y57i *= p2;
              y68r *= p2; y68i *= p2;
              y15r *= p3; y15i *= p3;
              y26r *= p3; y26i *= p3;
              y37r *= p3; y37i *= p3;
              y48r *= p3; y48i *= p3;
              gxr[j1] += y12r + y13r + y15r;
              gxr[j2] -= y12r - y24r - y26r;
              gxr[j3] += y34r - y13r + y37r;
              gxr[j4] -= y34r + y24r - y48r;
              gxr[j5] += y56r + y57r - y15r;
              gxr[j6] -= y56r - y68r + y26r;
              gxr[j7] += y78r - y57r - y37r;
              gxr[j8] -= y78r + y68r + y48r;
              gxi[j1] += y12i + y13i + y15i;
              gxi[j2] -= y12i - y24i - y26i;
              gxi[j3] += y34i - y13i + y37i;
              gxi[j4] -= y34i + y24i - y48i;
              gxi[j5] += y56i + y57i - y15i;
              gxi[j6] -= y56i - y68i + y26i;
              gxi[j7] += y78i - y57i - y37i;
              gxi[j8] -= y78i + y68i + y48i;
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            j8 = OFF3(0, i2,   i3);
            j6 = OFF3(0, i2-1, i3);
            j4 = OFF3(0, i2,   i3-1);
            j2 = OFF3(0, i2-1, i3-1);
            x2r = xr[j2]; x2i = xi[j2]; 
            x4r = xr[j4]; x4i = xi[j4];
            x6r = xr[j6]; x6i = xi[j6];
            x8r = xr[j8]; x8i = xi[j8];
            for (i1 = 1; i1 < n1; ++i1) {
              /* Peak the values at the 8 corners of the cube. */
              j1 = j2++; x1r = x2r; x2r = xr[j2]; x1i = x2i; x2i = xi[j2];
              j3 = j4++; x3r = x4r; x4r = xr[j4]; x3i = x4i; x4i = xi[j4];
              j5 = j6++; x5r = x6r; x6r = xr[j6]; x5i = x6i; x6i = xi[j6];
              j7 = j8++; x7r = x8r; x8r = xr[j8]; x7i = x8i; x8i = xi[j8];
              /* Compute the differences along all the 12 edges of the cube: */
              /*  - along 1st dim: */
              y12r = x1r - x2r; y12i = x1i - x2i;
              y34r = x3r - x4r; y34i = x3i - x4i;
              y56r = x5r - x6r; y56i = x5i - x6i;
              y78r = x7r - x8r; y78i = x7i - x8i;
              /*  - along 2nd dim: */
              y13r = x1r - x3r; y13i = x1i - x3i;
              y24r = x2r - x4r; y24i = x2i - x4i;
              y57r = x5r - x7r; y57i = x5i - x7i;
              y68r = x6r - x8r; y68i = x6i - x8i;
              /*  - along 3rd dim: */
              y15r = x1r - x5r; y15i = x1i - x5i;
              y26r = x2r - x6r; y26i = x2i - x6i;
              y37r = x3r - x7r; y37i = x3i - x7i;
              y48r = x4r - x8r; y48i = x4i - x8i;
              /* Compute the cost and integrate its gradient. */
              r = ((SQ(y12r) + SQ(y34r) + SQ(y56r) + SQ(y78r))*q1 +
                   (SQ(y13r) + SQ(y24r) + SQ(y57r) + SQ(y68r))*q2 +
                   (SQ(y15r) + SQ(y26r) + SQ(y37r) + SQ(y48r))*q3 +
                   (SQ(y12i) + SQ(y34i) + SQ(y56i) + SQ(y78i))*q1 +
                   (SQ(y13i) + SQ(y24i) + SQ(y57i) + SQ(y68i))*q2 +
                   (SQ(y15i) + SQ(y26i) + SQ(y37i) + SQ(y48i))*q3);
              fx += sqrt(r + s);
            }
          }
        }
      }
    }

  }

  /* Remove the "bias" and make sure the result is non-negative (it can only
     be negative due to rounding errors). */
  fx -= (n1 - 1)*(n2 - 1)*(n3 - 1)*eps;
  if (fx < 0.0) {
    fx = 0.0;
  }
  return fx;
}

/*
 * NOTATIONS FOR 4-D VOLUME X(i1,i2,i3,i4)
 *
 *                 i3  i2
 *                  | /
 *                  |/              X1 = X(i1-1,i2-1,i3-1,i4-1)
 *        X7--------X8---> i1       X2 = X(i1  ,i2-1,i3-1,i4-1)
 *       /:        /|               X3 = X(i1-1,i2  ,i3-1,i4-1)
 *      / :       / |               X4 = X(i1  ,i2  ,i3-1,i4-1)
 *     X5--------X6 |               X5 = X(i1-1,i2-1,i3  ,i4-1)
 *     |  X3.....|..X4              X6 = X(i1  ,i2-1,i3  ,i4-1)
 *     | '       | /                X7 = X(i1-1,i2  ,i3  ,i4-1)
 *     |'        |/                 X8 = X(i1  ,i2  ,i3  ,i4-1)
 *     X1--------X2
 *     :
 *     :               i3  i2
 *     :                | /
 *      :               |/              X9  = X(i1-1,i2-1,i3-1,i4  )
 *      :     X15-------X16---> i1      X10 = X(i1  ,i2-1,i3-1,i4  )
 *   i4 :     /:        /|              X11 = X(i1-1,i2  ,i3-1,i4  )
 *       :   / :       / |              X12 = X(i1  ,i2  ,i3-1,i4  )
 *       :  X13-------X14|              X13 = X(i1-1,i2-1,i3  ,i4  )
 *       :  |  X11....|..X12            X14 = X(i1  ,i2-1,i3  ,i4  )
 *        : | '       | /               X15 = X(i1-1,i2  ,i3  ,i4  )
 *        : |'        |/                X16 = X(i1  ,i2  ,i3  ,i4  )
 *        v X9--------X10
 *
 */
double rgl_tv4d(const double x[],
                const long n1, const long n2, const long n3, const long n4,
                const double w1, const double w2, const double w3, const double w4,
                const double eps, double gx[], unsigned int flags)
{
  double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16;
  double r, s, fx;
  double y12, y34, y56, y78, y910, y1112, y1314, y1516;
  double y13, y24, y57, y68, y911, y1012, y1315, y1416;
  double y15, y26, y37, y48, y913, y1014, y1115, y1216;
  double y19, y210, y311, y412, y513, y614, y715, y816;

  long i1, i2, i3, i4;
  long j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16;
  int compute_gradient;

  compute_gradient = ((flags & COMPUTE_GRADIENT) != 0);
  if (x == NULL || n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 ||
      w1 < 0.0 || w2 < 0.0 || w3 < 0.0 || w4 < 0.0 || eps < 0.0 ||
      (gx == NULL && compute_gradient)) {
    return -1.0;
  }
  if ((flags & RGL_STORE_GRADIENT) != 0) {
    memset(gx, 0, n1*n2*n3*n4*sizeof(gx[0]));
  }
  fx = 0.0;
  s = eps*eps;

  if (METHOD(flags, RGL_TOTVAR_FORWARD)) {
    /*
     * The squared norm of the spatial gradient is the sum of the squared
     * forward differences along each direction (standard discretization).
     */

    if (w1 == w2 && w2 == w3 && w3 == w4) /* same weights along all directions */ {
      double p, q = w1;
      if (compute_gradient) {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j3 = OFF4(0, i2,   i3-1, i4-1);
              j5 = OFF4(0, i2-1, i3  , i4-1);
              j9 = OFF4(0, i2-1, i3-1, i4  );
              x2 = x[j2];
              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 4 corners of the cube */
                j1 = j2++;
                x1 = x2;
                x2 = x[j2];
                x3 = x[j3];
                x5 = x[j5];
                x9 = x[j9];
                /* Compute the differences along the 3 directions: */
                y12 = x1 - x2;
                y13 = x1 - x3;
                y15 = x1 - x5;
                y19 = x1 - x9;
                /* Compute the cost and integrate its gradient. */
                r = (SQ(y12) + SQ(y13) + SQ(y15) + SQ(y19))*q;
                p = sqrt(r + s);
                fx += p;
                p = q/p;
                gx[j1] += (y12 + y13 + y15 + y19)*p;
                gx[j2] -= y12*p;
                gx[j3] -= y13*p;
                gx[j5] -= y15*p;
                gx[j9] -= y19*p;
                ++j3;
                ++j5;
                ++j9;
              }
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j3 = OFF4(0, i2,   i3-1, i4-1);
              j5 = OFF4(0, i2-1, i3  , i4-1);
              j9 = OFF4(0, i2-1, i3-1, i4  );
              x2 = x[j2];
              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 4 corners of the cube */
                x1 = x2;
                x2 = x[++j2];
                x3 = x[j3++];
                x5 = x[j5++];
                x9 = x[j9++];
                /* Compute the differences along the 3 directions: */
                y12 = x1 - x2;
                y13 = x1 - x3;
                y15 = x1 - x5;
                y19 = x1 - x9;
                /* Compute the cost and integrate its gradient. */
                r = (SQ(y12) + SQ(y13) + SQ(y15) + SQ(y19))*q;
                fx += sqrt(r + s);
              }
            }
          }
        }
      }
    } else /* not same weights along all directions */ {
      double q1 = w1;
      double q2 = w2;
      double q3 = w3;
      double q4 = w4;
      if (compute_gradient) {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j3 = OFF4(0, i2,   i3-1, i4-1);
              j5 = OFF4(0, i2-1, i3  , i4-1);
              j9 = OFF4(0, i2-1, i3-1, i4  );
              x2 = x[j2];
              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 4 corners of the cube */
                j1 = j2++;
                x1 = x2;
                x2 = x[j2];
                x3 = x[j3];
                x5 = x[j5];
                x9 = x[j9];
                /* Compute the differences along the 3 directions: */
                y12 = x1 - x2;
                y13 = x1 - x3;
                y15 = x1 - x5;
                y19 = x1 - x9;
                /* Compute the cost and integrate its gradient. */
                r =  SQ(y12)*q1 + SQ(y13)*q2 + SQ(y15)*q3 + SQ(y19)*q4;
                r = sqrt(r + s);
                fx += r;
                r = 1.0/r;
                y12 *= r*q1;
                y13 *= r*q2;
                y15 *= r*q3;
                y19 *= r*q4;
                gx[j1] += y12 + y13 + y15 + y19;
                gx[j2] -= y12;
                gx[j3] -= y13;
                gx[j5] -= y15;
                gx[j9] -= y19;
                ++j3;
                ++j5;
                ++j9;
              }
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j3 = OFF4(0, i2,   i3-1, i4-1);
              j5 = OFF4(0, i2-1, i3  , i4-1);
              j9 = OFF4(0, i2-1, i3-1, i4  );
              x2 = x[j2];
              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at the 8 corners of the cube */
                x1 = x2;
                x2 = x[++j2];
                x3 = x[j3++];
                x5 = x[j5++];
                x9 = x[j9++];
                /* Compute the differences along the 3 directions: */
                y12 = x1 - x2;
                y13 = x1 - x3;
                y15 = x1 - x5;
                y19 = x1 - x9;
                /* Compute the cost and integrate its gradient. */
                r = SQ(y12)*q1 + SQ(y13)*q2 + SQ(y15)*q3 + SQ(y19)*q4;
                fx += sqrt(r + s);
              }
            }
          }
        }
      }
    }
  } else /* same as RGL_TOTVAR_ISOTROPIC */ {
    if (w1 == w2 && w2 == w3 && w3 == w4) /* same weights along all directions */ {
      double p, q = w1/8.0;
      if (compute_gradient) {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j8 = OFF4(0, i2,   i3,   i4-1);
              j6 = OFF4(0, i2-1, i3,   i4-1);
              j4 = OFF4(0, i2,   i3-1, i4-1);
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j16 = OFF4(0, i2,   i3,   i4);
              j14 = OFF4(0, i2-1, i3,   i4);
              j12 = OFF4(0, i2,   i3-1, i4);
              j10 = OFF4(0, i2-1, i3-1, i4);

              x2 = x[j2];
              x4 = x[j4];
              x6 = x[j6];
              x8 = x[j8];
              x10 = x[j10];
              x12 = x[j12];
              x14 = x[j14];
              x16 = x[j16];

              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 8 corners of the hypercube */
                j1 = j2++;   x1 = x2;   x2 = x[j2];
                j3 = j4++;   x3 = x4;   x4 = x[j4];
                j5 = j6++;   x5 = x6;   x6 = x[j6];
                j7 = j8++;   x7 = x8;   x8 = x[j8];
                j9 = j10++;  x9 = x10;  x10 = x[j10];
                j11 = j12++; x11 = x12; x12 = x[j12];
                j13 = j14++; x13 = x14; x14 = x[j14];
                j15 = j16++; x15 = x16; x16 = x[j16];

                /* Compute the differences along all the 24 edges of the hypercube: */
                /*  - along 1st dim: */
                y12 = x1 - x2;
                y34 = x3 - x4;
                y56 = x5 - x6;
                y78 = x7 - x8;
                y910 = x9 - x10;
                y1112 = x11 - x12;
                y1314 = x13 - x14;
                y1516 = x15 - x16;
                /*  - along 2nd dim: */
                y13 = x1 - x3;
                y24 = x2 - x4;
                y57 = x5 - x7;
                y68 = x6 - x8;
                y911 = x9 - x11;
                y1012 = x10 - x12;
                y1315 = x13 - x15;
                y1416 = x14 - x16;
                /*  - along 3rd dim: */
                y15 = x1 - x5;
                y26 = x2 - x6;
                y37 = x3 - x7;
                y48 = x4 - x8;
                y913 = x9 - x13;
                y1014 = x10 - x14;
                y1115 = x11 - x15;
                y1216 = x12 - x16;
                /*  - along 4th dim: */
                y19 = x1 - x9;
                y210 = x2 - x10;
                y311 = x3 - x11;
                y412 = x4 - x12;
                y513 = x5 - x13;
                y614 = x6 - x14;
                y715 = x7 - x15;
                y816 = x8 - x16;

                /* Compute the cost and integrate its gradient. */
                r = (SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78) + SQ(y910) + SQ(y1112) + SQ(y1314) + SQ(y1516) +
                     SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68) + SQ(y911) + SQ(y1012) + SQ(y1315) + SQ(y1416) +
                     SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48) + SQ(y913) + SQ(y1014) + SQ(y1115) + SQ(y1216) +
                     SQ(y19) + SQ(y210) + SQ(y311) + SQ(y412) + SQ(y513) + SQ(y614) + SQ(y715) + SQ(y816))*q;
                p = sqrt(r + s);
                fx += p;
                p = q/p;
                gx[j1] += (y12 + y13 + y15 + y19)*p;
                gx[j2] -= (y12 - y24 - y26 - y210)*p;
                gx[j3] += (y34 - y13 + y37 + y311)*p;
                gx[j4] -= (y34 + y24 - y48 - y412)*p;
                gx[j5] += (y56 + y57 - y15 + y513)*p;
                gx[j6] -= (y56 - y68 + y26 - y614)*p;
                gx[j7] += (y78 - y57 - y37 + y715)*p;
                gx[j8] -= (y78 + y68 + y48 - y816)*p;
                gx[j9] += (y910 + y911 + y913 - y19)*p;
                gx[j10] -= (y910 - y1012 - y1014 + y210)*p;
                gx[j11] += (y1112 - y911 + y1115 - y311)*p;
                gx[j12] -= (y1112 + y1012 - y1216 + y412)*p;
                gx[j13] += (y1314 + y1315 - y913 - y513)*p;
                gx[j14] -= (y1314 - y1416 + y1014 + y614)*p;
                gx[j15] += (y1516 - y1315 - y1115 - y715)*p;
                gx[j16] -= (y1516 + y1416 + y1216 + y816)*p;
              }
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j8 = OFF4(0, i2,   i3,   i4-1);
              j6 = OFF4(0, i2-1, i3,   i4-1);
              j4 = OFF4(0, i2,   i3-1, i4-1);
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j16 = OFF4(0, i2,   i3,   i4);
              j14 = OFF4(0, i2-1, i3,   i4);
              j12 = OFF4(0, i2,   i3-1, i4);
              j10 = OFF4(0, i2-1, i3-1, i4);

              x2 = x[j2];
              x4 = x[j4];
              x6 = x[j6];
              x8 = x[j8];
              x10 = x[j10];
              x12 = x[j12];
              x14 = x[j14];
              x16 = x[j16];

              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 8 corners of the hypercube */
                j1 = j2++;   x1 = x2;   x2 = x[j2];
                j3 = j4++;   x3 = x4;   x4 = x[j4];
                j5 = j6++;   x5 = x6;   x6 = x[j6];
                j7 = j8++;   x7 = x8;   x8 = x[j8];
                j9 = j10++;  x9 = x10;  x10 = x[j10];
                j11 = j12++; x11 = x12; x12 = x[j12];
                j13 = j14++; x13 = x14; x14 = x[j14];
                j15 = j16++; x15 = x16; x16 = x[j16];

                /* Compute the differences along all the 24 edges of the hypercube: */
                /*  - along 1st dim: */
                y12 = x1 - x2;
                y34 = x3 - x4;
                y56 = x5 - x6;
                y78 = x7 - x8;
                y910 = x9 - x10;
                y1112 = x11 - x12;
                y1314 = x13 - x14;
                y1516 = x15 - x16;
                /*  - along 2nd dim: */
                y13 = x1 - x3;
                y24 = x2 - x4;
                y57 = x5 - x7;
                y68 = x6 - x8;
                y911 = x9 - x11;
                y1012 = x10 - x12;
                y1315 = x13 - x15;
                y1416 = x14 - x16;
                /*  - along 3rd dim: */
                y15 = x1 - x5;
                y26 = x2 - x6;
                y37 = x3 - x7;
                y48 = x4 - x8;
                y913 = x9 - x13;
                y1014 = x10 - x14;
                y1115 = x11 - x15;
                y1216 = x12 - x16;
                /*  - along 4th dim: */
                y19 = x1 - x9;
                y210 = x2 - x10;
                y311 = x3 - x11;
                y412 = x4 - x12;
                y513 = x5 - x13;
                y614 = x6 - x14;
                y715 = x7 - x15;
                y816 = x8 - x16;

                /* Compute the cost and integrate its gradient. */
                r = (SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78) + SQ(y910) + SQ(y1112) + SQ(y1314) + SQ(y1516) +
                     SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68) + SQ(y911) + SQ(y1012) + SQ(y1315) + SQ(y1416) +
                     SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48) + SQ(y913) + SQ(y1014) + SQ(y1115) + SQ(y1216) +
                     SQ(y19) + SQ(y210) + SQ(y311) + SQ(y412) + SQ(y513) + SQ(y614) + SQ(y715) + SQ(y816))*q;
                fx += sqrt(r + s);
              }
            }
          }
        }
      }
    } else /* not same weights along all directions */ {
      double p1, q1 = w1/8.0;
      double p2, q2 = w2/8.0;
      double p3, q3 = w3/8.0;
      double p4, q4 = w4/8.0;
      if (compute_gradient) {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j8 = OFF4(0, i2,   i3,   i4-1);
              j6 = OFF4(0, i2-1, i3,   i4-1);
              j4 = OFF4(0, i2,   i3-1, i4-1);
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j16 = OFF4(0, i2,   i3,   i4);
              j14 = OFF4(0, i2-1, i3,   i4);
              j12 = OFF4(0, i2,   i3-1, i4);
              j10 = OFF4(0, i2-1, i3-1, i4);

              x2 = x[j2];
              x4 = x[j4];
              x6 = x[j6];
              x8 = x[j8];
              x10 = x[j10];
              x12 = x[j12];
              x14 = x[j14];
              x16 = x[j16];

              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 8 corners of the hypercube */
                j1 = j2++;   x1 = x2;   x2 = x[j2];
                j3 = j4++;   x3 = x4;   x4 = x[j4];
                j5 = j6++;   x5 = x6;   x6 = x[j6];
                j7 = j8++;   x7 = x8;   x8 = x[j8];
                j9 = j10++;  x9 = x10;  x10 = x[j10];
                j11 = j12++; x11 = x12; x12 = x[j12];
                j13 = j14++; x13 = x14; x14 = x[j14];
                j15 = j16++; x15 = x16; x16 = x[j16];

                /* Compute the differences along all the 24 edges of the hypercube: */
                /*  - along 1st dim: */
                y12 = x1 - x2;
                y34 = x3 - x4;
                y56 = x5 - x6;
                y78 = x7 - x8;
                y910 = x9 - x10;
                y1112 = x11 - x12;
                y1314 = x13 - x14;
                y1516 = x15 - x16;
                /*  - along 2nd dim: */
                y13 = x1 - x3;
                y24 = x2 - x4;
                y57 = x5 - x7;
                y68 = x6 - x8;
                y911 = x9 - x11;
                y1012 = x10 - x12;
                y1315 = x13 - x15;
                y1416 = x14 - x16;
                /*  - along 3rd dim: */
                y15 = x1 - x5;
                y26 = x2 - x6;
                y37 = x3 - x7;
                y48 = x4 - x8;
                y913 = x9 - x13;
                y1014 = x10 - x14;
                y1115 = x11 - x15;
                y1216 = x12 - x16;
                /*  - along 4th dim: */
                y19 = x1 - x9;
                y210 = x2 - x10;
                y311 = x3 - x11;
                y412 = x4 - x12;
                y513 = x5 - x13;
                y614 = x6 - x14;
                y715 = x7 - x15;
                y816 = x8 - x16;

                /* Compute the cost and integrate its gradient. */
                r = ((SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78) + SQ(y910) + SQ(y1112) + SQ(y1314) + SQ(y1516))*q1 +
                     (SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68) + SQ(y911) + SQ(y1012) + SQ(y1315) + SQ(y1416))*q2 +
                     (SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48) + SQ(y913) + SQ(y1014) + SQ(y1115) + SQ(y1216))*q3 +
                     (SQ(y19) + SQ(y210) + SQ(y311) + SQ(y412) + SQ(y513) + SQ(y614) + SQ(y715) + SQ(y816))*q4);
                r = sqrt(r + s);
                fx += r;
                r = 1.0/r;
                p1 = r*q1;
                p2 = r*q2;
                p3 = r*q3;
                p4 = r*q4;
                y12 *= p1;
                y34 *= p1;
                y56 *= p1;
                y78 *= p1;
                y910 *= p1;
                y1112 *= p1;
                y1314 *= p1;
                y1516 *= p1;
                y13 *= p2;
                y24 *= p2;
                y57 *= p2;
                y68 *= p2;
                y911 *= p2;
                y1012 *= p2;
                y1315 *= p2;
                y1416 *= p2;
                y15 *= p3;
                y26 *= p3;
                y37 *= p3;
                y48 *= p3;
                y913 *= p3;
                y1014 *= p3;
                y1115 *= p3;
                y1216 *= p3;
                y19 *= p4;
                y210 *= p4;
                y311 *= p4;
                y412 *= p4;
                y513 *= p4;
                y614 *= p4;
                y715 *= p4;
                y816 *= p4;
                gx[j1] += y12 + y13 + y15 + y19;
                gx[j2] -= y12 - y24 - y26 - y210;
                gx[j3] += y34 - y13 + y37 + y311;
                gx[j4] -= y34 + y24 - y48 - y412;
                gx[j5] += y56 + y57 - y15 + y513;
                gx[j6] -= y56 - y68 + y26 - y614;
                gx[j7] += y78 - y57 - y37 + y715;
                gx[j8] -= y78 + y68 + y48 - y816;
                gx[j9] += y910 + y911 + y913 - y19;
                gx[j10] -= y910 - y1012 - y1014 + y210;
                gx[j11] += y1112 - y911 + y1115 - y311;
                gx[j12] -= y1112 + y1012 - y1216 + y412;
                gx[j13] += y1314 + y1315 - y913 - y513;
                gx[j14] -= y1314 - y1416 + y1014 + y614;
                gx[j15] += y1516 - y1315 - y1115 - y715;
                gx[j16] -= y1516 + y1416 + y1216 + y816;
              }
            }
          }
        }
      } else /* do not compute gradients */ {
        for (i4 = 1; i4 < n4; ++i4) {
          for (i3 = 1; i3 < n3; ++i3) {
            for (i2 = 1; i2 < n2; ++i2) {
              j8 = OFF4(0, i2,   i3,   i4-1);
              j6 = OFF4(0, i2-1, i3,   i4-1);
              j4 = OFF4(0, i2,   i3-1, i4-1);
              j2 = OFF4(0, i2-1, i3-1, i4-1);
              j16 = OFF4(0, i2,   i3,   i4);
              j14 = OFF4(0, i2-1, i3,   i4);
              j12 = OFF4(0, i2,   i3-1, i4);
              j10 = OFF4(0, i2-1, i3-1, i4);

              x2 = x[j2];
              x4 = x[j4];
              x6 = x[j6];
              x8 = x[j8];
              x10 = x[j10];
              x12 = x[j12];
              x14 = x[j14];
              x16 = x[j16];

              for (i1 = 1; i1 < n1; ++i1) {
                /* Peak the values at 8 corners of the hypercube */
                j1 = j2++;   x1 = x2;   x2 = x[j2];
                j3 = j4++;   x3 = x4;   x4 = x[j4];
                j5 = j6++;   x5 = x6;   x6 = x[j6];
                j7 = j8++;   x7 = x8;   x8 = x[j8];
                j9 = j10++;  x9 = x10;  x10 = x[j10];
                j11 = j12++; x11 = x12; x12 = x[j12];
                j13 = j14++; x13 = x14; x14 = x[j14];
                j15 = j16++; x15 = x16; x16 = x[j16];

                /* Compute the differences along all the 24 edges of the hypercube: */
                /*  - along 1st dim: */
                y12 = x1 - x2;
                y34 = x3 - x4;
                y56 = x5 - x6;
                y78 = x7 - x8;
                y910 = x9 - x10;
                y1112 = x11 - x12;
                y1314 = x13 - x14;
                y1516 = x15 - x16;
                /*  - along 2nd dim: */
                y13 = x1 - x3;
                y24 = x2 - x4;
                y57 = x5 - x7;
                y68 = x6 - x8;
                y911 = x9 - x11;
                y1012 = x10 - x12;
                y1315 = x13 - x15;
                y1416 = x14 - x16;
                /*  - along 3rd dim: */
                y15 = x1 - x5;
                y26 = x2 - x6;
                y37 = x3 - x7;
                y48 = x4 - x8;
                y913 = x9 - x13;
                y1014 = x10 - x14;
                y1115 = x11 - x15;
                y1216 = x12 - x16;
                /*  - along 4th dim: */
                y19 = x1 - x9;
                y210 = x2 - x10;
                y311 = x3 - x11;
                y412 = x4 - x12;
                y513 = x5 - x13;
                y614 = x6 - x14;
                y715 = x7 - x15;
                y816 = x8 - x16;

                /* Compute the cost and integrate its gradient. */
                r = ((SQ(y12) + SQ(y34) + SQ(y56) + SQ(y78) + SQ(y910) + SQ(y1112) + SQ(y1314) + SQ(y1516))*q1 +
                     (SQ(y13) + SQ(y24) + SQ(y57) + SQ(y68) + SQ(y911) + SQ(y1012) + SQ(y1315) + SQ(y1416))*q2 +
                     (SQ(y15) + SQ(y26) + SQ(y37) + SQ(y48) + SQ(y913) + SQ(y1014) + SQ(y1115) + SQ(y1216))*q3 +
                     (SQ(y19) + SQ(y210) + SQ(y311) + SQ(y412) + SQ(y513) + SQ(y614) + SQ(y715) + SQ(y816))*q4);
                fx += sqrt(r + s);
              }
            }
          }
        }
      }
    }
  }

  /* Remove the "bias" and make sure the result is non-negative (it can only
     be negative due to rounding errors). */
  fx -= (n1 - 1)*(n2 - 1)*(n3 - 1)*(n4 - 1)*eps;
  if (fx < 0.0) {
    fx = 0.0;
  }
  return fx;
}
