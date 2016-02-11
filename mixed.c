/*
 * mixed.c --
 *
 * Implements N-D + t edge-preserving regularization (for N = 2 or 3).
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

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Inline functions.  */
static inline double square(const double) __attribute__((always_inline));
static inline double square(const double x) { return x*x; }

#define INDEX3(i1,i2,i3)    ((i1) + n1*((i2) + n2*(i3)))
#define INDEX4(i1,i2,i3,i4) ((i1) + n1*((i2) + n2*((i3) + n3*(i4))))
#define X(a)                x[j##a]
#define G(a)                g[j##a]
#define SQDF(a,b)           square(X(a) - X(b)) /* squared difference */
#define GRD2(a,b,c,r)       G(a) += (2*X(a) - X(b) - X(c))*(r)
#define GRD3(a,b,c,d,r)     G(a) += (3*X(a) - X(b) - X(c) - X(d))*(r)

#ifndef M_SQRT2
static const double M_SQRT2 = 1.414213562373095048801688724209698079;
#endif

double mixed_regul_2dpt(double mu1, double eps1,
                        double mu2, double eps2,
                        long n1, long n2, long n3,
                        const double* x, double *g,
                        int clr)
{
  double s1, s2, f1, f2;
  long i1, i2, i3;

  /* Account for the increased number of finite differences to compute the
     squared norm of the gradient. */
  mu1  /= M_SQRT2;
  eps1 *= M_SQRT2;

  /* Initialize and loop over elements. */
  s1 = eps1*eps1;
  s2 = eps2*eps2;
  f1 = 0;
  f2 = 0;
  if (g != NULL) {
    if (clr) {
      memset(g, 0, n1*n2*n3*sizeof(g[0]));
    }
    for (i3 = 0; i3 < n3; ++i3) {
      for (i2 = 1; i2 < n2; ++i2) {
        long j01 = INDEX3(0,i2-1,i3);
        long j11 = INDEX3(0,i2  ,i3);
        for (i1 = 1; i1 < n1; ++i1) {
          long j00 = j01++;
          long j10 = j11++;
          double r = sqrt(SQDF(00,01) +
                          SQDF(10,11) +
                          SQDF(00,10) +
                          SQDF(01,11) + s1);
          f1 += r;
          r = mu1/r;
          GRD2(00,01,10,r);
          GRD2(01,00,11,r);
          GRD2(10,00,11,r);
          GRD2(11,01,10,r);
        }
      }
      if (i3 >= 1) {
        long j0 = INDEX3(0,0,i3-1);
        long j1 = INDEX3(0,0,i3);
        for (i2 = 0; i2 < n2; ++i2) {
          for (i1 = 0; i1 < n1; ++i1, ++j0, ++j1) {
            double t = x[j1] - x[j0];
            double r = sqrt(t*t + s2);
            f2 += r;
            t *= mu2/r;
            g[j0] -= t;
            g[j1] += t;
          }
        }
      }
    }
  } else {
    /* no gradient */
    for (i3 = 0; i3 < n3; ++i3) {
      for (i2 = 1; i2 < n2; ++i2) {
        long j01 = INDEX3(0,i2-1,i3);
        long j11 = INDEX3(0,i2  ,i3);
        for (i1 = 1; i1 < n1; ++i1) {
          long j00 = j01++;
          long j10 = j11++;
          f1 += sqrt(SQDF(00,01) +
                     SQDF(10,11) +
                     SQDF(00,10) +
                     SQDF(01,11) + s1);
        }
      }
      if (i3 >= 1) {
        long j0 = INDEX3(0,0,i3-1);
        long j1 = INDEX3(0,0,i3);
        for (i2 = 0; i2 < n2; ++i2) {
          for (i1 = 0; i1 < n1; ++i1, ++j0, ++j1) {
            double t = x[j1] - x[j0];
            f2 += sqrt(t*t + s2);
          }
        }
      }
    }
  }
  return (mu1*(f1 - (n1 - 1)*(n2 - 1)*n3*eps1) +
          mu2*(f2 - n1*n2*(n3 - 1)*eps2));
}

double mixed_regul_3dpt(double mu1, double eps1,
                        double mu2, double eps2,
                        long n1, long n2, long n3, long n4,
                        const double* x, double *g,
                        int clr)
{
  double s1, s2;
  double f1, f2 = 0;
  long i1, i2, i3, i4;

  /* Account for the increased number of finite differences to compute the
     squared norm of the gradient. */
  mu1  /= 2;
  eps1 *= 2;

  /* Initialize and loop over elements. */
  s1 = eps1*eps1;
  s2 = eps2*eps2;
  f1 = 0;
  f2 = 0;
  if (g != NULL) {
    if (clr) {
      memset(g, 0, n1*n2*n3*n4*sizeof(g[0]));
    }
    for (i4 = 0; i4 < n4; ++i4) {
      for (i3 = 1; i3 < n3; ++i3) {
        for (i2 = 1; i2 < n2; ++i2) {
          long j001 = INDEX4(0,i2-1,i3-1,i4);
          long j011 = INDEX4(0,i2  ,i3-1,i4);
          long j101 = INDEX4(0,i2-1,i3  ,i4);
          long j111 = INDEX4(0,i2  ,i3  ,i4);
          for (i1 = 1; i1 < n1; ++i1) {
            long j000 = j001++;
            long j010 = j011++;
            long j100 = j101++;
            long j110 = j111++;
            double r = sqrt(SQDF(001,000) +
                            SQDF(011,010) +
                            SQDF(101,100) +
                            SQDF(111,110) +
                            SQDF(010,000) +
                            SQDF(011,001) +
                            SQDF(110,100) +
                            SQDF(111,101) +
                            SQDF(100,000) +
                            SQDF(101,001) +
                            SQDF(110,010) +
                            SQDF(111,011) + s1);
            f1 += r;
            r = mu1/r;
            GRD3(000,001,010,100,r);
            GRD3(001,000,011,101,r);
            GRD3(010,011,000,110,r);
            GRD3(011,010,001,111,r);
            GRD3(100,101,110,000,r);
            GRD3(101,100,111,001,r);
            GRD3(110,111,100,010,r);
            GRD3(111,110,101,011,r);
          }
        }
      }
      if (i4 >= 1) {
        long j0 = INDEX4(0,0,0,i4-1);
        long j1 = INDEX4(0,0,0,i4);
        for (i3 = 0; i3 < n3; ++i3) {
          for (i2 = 0; i2 < n2; ++i2) {
            for (i1 = 0; i1 < n1; ++i1, ++j0, ++j1) {
              double t = x[j1] - x[j0];
              double r = sqrt(t*t + s2);
              f2 += r;
              t *= mu2/r;
              g[j0] -= t;
              g[j1] += t;
            }
          }
        }
      }
    }
  } else {
    /* no gradient */
    for (i4 = 0; i4 < n4; ++i4) {
      for (i3 = 1; i3 < n3; ++i3) {
        for (i2 = 1; i2 < n2; ++i2) {
          long j001 = INDEX4(0,i2-1,i3-1,i4);
          long j011 = INDEX4(0,i2  ,i3-1,i4);
          long j101 = INDEX4(0,i2-1,i3  ,i4);
          long j111 = INDEX4(0,i2  ,i3  ,i4);
          for (i1 = 1; i1 < n1; ++i1) {
            long j000 = j001++;
            long j010 = j011++;
            long j100 = j101++;
            long j110 = j111++;
            f1 += sqrt(SQDF(001,000) +
                       SQDF(011,010) +
                       SQDF(101,100) +
                       SQDF(111,110) +
                       SQDF(010,000) +
                       SQDF(011,001) +
                       SQDF(110,100) +
                       SQDF(111,101) +
                       SQDF(100,000) +
                       SQDF(101,001) +
                       SQDF(110,010) +
                       SQDF(111,011) + s1);
          }
        }
      }
      if (i4 >= 1) {
        long j0 = INDEX4(0,0,0,i4-1);
        long j1 = INDEX4(0,0,0,i4);
        for (i3 = 0; i3 < n3; ++i3) {
          for (i2 = 0; i2 < n2; ++i2) {
            for (i1 = 0; i1 < n1; ++i1, ++j0, ++j1) {
              double t = x[j1] - x[j0];
              f2 += sqrt(t*t + s2);
            }
          }
        }
      }
    }
  }
  return (mu1*(f1 - (n1 - 1)*(n2 - 1)*(n3 - 1)*n4*eps1) +
          mu2*(f2 - n1*n2*n3*(n4 - 1)*eps2));
}
