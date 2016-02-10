
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define INDEX3(i1,i2,i3)    ((i1) + n1*((i2) + n2*(i3)))
#define INDEX4(i1,i2,i3,i4) ((i1) + n1*((i2) + n2*((i3) + n3*(i4))))
/*
#define X(i1,i2,i3,i4)      x[INDEX4(i1,i2,i3,i4)]
#define G(i1,i2,i3,i4)      g[INDEX4(i1,i2,i3,i4)]
*/
#define SQUARE(a)           ((a)*(a))
#define X(a)                x[j##a]
#define G(a)                g[j##a]
#define Q(a,b)              SQUARE(x##a##_x##b)
#define DIF(a,b)            double x##a##_x##b = X(a) - X(b)
#define GRD2(a,b,c)         G(a) += (2*X(a) - X(b) - X(c))*r
#define GRD3(a,b,c,d)       G(a) += (3*X(a) - X(b) - X(c) - X(d))*r

double mixed_regul_2dpt(double mu1, double eps1,
                        double mu2, double eps2,
                        long n1, long n2, long n3,
                        const double* x, double *g,
                        int clr)
{
  double s1 = eps1*eps1;
  double s2 = eps2*eps2;
  double r, f1 = 0, f2 = 0;
  long i1, i2, i3;

  if (g != NULL) {
    if (clr && g != NULL) {
      memset(g, 0, n1*n2*n3*sizeof(g[0]));
    }
    for (i3 = 0; i3 < n3; ++i3) {
      if (i3 >= 1) {
        for (i2 = 1; i2 < n2; ++i2) {
          for (i1 = 1; i1 < n1; ++i1) {
            long j1 = INDEX3(i1,i2,i3);
            long j0 = j1 - n2; /* INDEX3(i1,i2,i3-1) */
            double x1_x0 = x[j1] - x[j0];
            r = sqrt(SQUARE(x1_x0) + s2);
            f2 += r;
            r = mu2/r;
            g[j0] -= r*x1_x0;
            g[j1] += r*x1_x0;
          }
        }
      }
      for (i2 = 1; i2 < n2; ++i2) {
        for (i1 = 1; i1 < n1; ++i1) {
          long j00 = INDEX3(i1-1,i2-1,i3);
          long j01 = INDEX3(i1  ,i2-1,i3);
          long j10 = INDEX3(i1-1,i2  ,i3);
          long j11 = INDEX3(i1  ,i2  ,i3);
          DIF(01,00);
          DIF(11,10);
          DIF(10,00);
          DIF(11,01);
          r = sqrt(Q(01,00) +
                   Q(11,10) +
                   Q(10,00) +
                   Q(11,01) + s1);
          f1 += r;
          r = mu1/r;
          GRD2(00,01,10);
          GRD2(01,00,11);
          GRD2(10,00,11);
          GRD2(11,01,10);
        }
      }
    }
  } else {
    /* no gradient */
    for (i3 = 0; i3 < n3; ++i3) {
      if (i3 >= 1) {
        for (i2 = 1; i2 < n2; ++i2) {
          for (i1 = 1; i1 < n1; ++i1) {
            long j1 = INDEX3(i1,i2,i3);
            long j0 = j1 - n2; /* INDEX3(i1,i2,i3-1) */
            double x1_x0 = x[j1] - x[j0];
            f2 += sqrt(SQUARE(x1_x0) + s2);
          }
        }
      }
      for (i2 = 1; i2 < n2; ++i2) {
        for (i1 = 1; i1 < n1; ++i1) {
          long j00 = INDEX3(i1-1,i2-1,i3);
          long j01 = INDEX3(i1  ,i2-1,i3);
          long j10 = INDEX3(i1-1,i2  ,i3);
          long j11 = INDEX3(i1  ,i2  ,i3);
          DIF(01,00);
          DIF(11,10);
          DIF(10,00);
          DIF(11,01);
          f1 += sqrt(Q(01,00) +
                     Q(11,10) +
                     Q(10,00) +
                     Q(11,01) + s1);
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
  double s1 = eps1*eps1;
  double s2 = eps2*eps2;
  double r, f1 = 0, f2 = 0;
  long i1, i2, i3, i4;

  if (g != NULL) {
    if (clr && g != NULL) {
      memset(g, 0, n1*n2*n3*n4*sizeof(g[0]));
    }
    for (i4 = 0; i4 < n4; ++i4) {
      if (i4 >= 1) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            for (i1 = 1; i1 < n1; ++i1) {
              long j1 = INDEX4(i1,i2,i3,i4);
              long j0 = j1 - n3; /* INDEX4(i1,i2,i3,i4-1) */
              double x1_x0 = x[j1] - x[j0];
              r = sqrt(SQUARE(x1_x0) + s2);
              f2 += r;
              r = mu2/r;
              g[j0] -= r*x1_x0;
              g[j1] += r*x1_x0;
            }
          }
        }
      }
      for (i3 = 1; i3 < n3; ++i3) {
        for (i2 = 1; i2 < n2; ++i2) {
          for (i1 = 1; i1 < n1; ++i1) {
            long j000 = INDEX4(i1-1,i2-1,i3-1,i4);
            long j001 = INDEX4(i1  ,i2-1,i3-1,i4);
            long j010 = INDEX4(i1-1,i2  ,i3-1,i4);
            long j011 = INDEX4(i1  ,i2  ,i3-1,i4);
            long j100 = INDEX4(i1-1,i2-1,i3  ,i4);
            long j101 = INDEX4(i1  ,i2-1,i3  ,i4);
            long j110 = INDEX4(i1-1,i2  ,i3  ,i4);
            long j111 = INDEX4(i1  ,i2  ,i3  ,i4);
            DIF(001,000);
            DIF(011,010);
            DIF(101,100);
            DIF(111,110);
            DIF(010,000);
            DIF(011,001);
            DIF(110,100);
            DIF(111,101);
            DIF(100,000);
            DIF(101,001);
            DIF(110,010);
            DIF(111,011);
            r = sqrt(Q(001,000) +
                     Q(011,010) +
                     Q(101,100) +
                     Q(111,110) +
                     Q(010,000) +
                     Q(011,001) +
                     Q(110,100) +
                     Q(111,101) +
                     Q(100,000) +
                     Q(101,001) +
                     Q(110,010) +
                     Q(111,011) + s1);
            f1 += r;
            r = mu1/r;
            GRD3(000,001,010,100);
            GRD3(001,000,011,101);
            GRD3(010,011,000,110);
            GRD3(011,010,001,111);
            GRD3(100,101,110,000);
            GRD3(101,100,111,001);
            GRD3(110,111,100,010);
            GRD3(111,110,101,011);
          }
        }
      }
    }
  } else {
    /* no gradient */
    for (i4 = 0; i4 < n4; ++i4) {
      if (i4 >= 1) {
        for (i3 = 1; i3 < n3; ++i3) {
          for (i2 = 1; i2 < n2; ++i2) {
            for (i1 = 1; i1 < n1; ++i1) {
              long j1 = INDEX4(i1,i2,i3,i4);
              long j0 = j1 - n3; /* INDEX4(i1,i2,i3,i4-1) */
              double x1_x0 = x[j1] - x[j0];
              f2 += sqrt(SQUARE(x1_x0) + s2);
            }
          }
        }
      }
      for (i3 = 1; i3 < n3; ++i3) {
        for (i2 = 1; i2 < n2; ++i2) {
          for (i1 = 1; i1 < n1; ++i1) {
            long j000 = INDEX4(i1-1,i2-1,i3-1,i4);
            long j001 = INDEX4(i1  ,i2-1,i3-1,i4);
            long j010 = INDEX4(i1-1,i2  ,i3-1,i4);
            long j011 = INDEX4(i1  ,i2  ,i3-1,i4);
            long j100 = INDEX4(i1-1,i2-1,i3  ,i4);
            long j101 = INDEX4(i1  ,i2-1,i3  ,i4);
            long j110 = INDEX4(i1-1,i2  ,i3  ,i4);
            long j111 = INDEX4(i1  ,i2  ,i3  ,i4);
            DIF(001,000);
            DIF(011,010);
            DIF(101,100);
            DIF(111,110);
            DIF(010,000);
            DIF(011,001);
            DIF(110,100);
            DIF(111,101);
            DIF(100,000);
            DIF(101,001);
            DIF(110,010);
            DIF(111,011);
            f1 += sqrt((Q(001,000) +
                        Q(011,010) +
                        Q(101,100) +
                        Q(111,110) +
                        Q(010,000) +
                        Q(011,001) +
                        Q(110,100) +
                        Q(111,101) +
                        Q(100,000) +
                        Q(101,001) +
                        Q(110,010) +
                        Q(111,011)) + s1);
          }
        }
      }
    }
  }
  return (mu1*(f1 - (n1 - 1)*(n2 - 1)*(n3 - 1)*n4*eps1) +
          mu2*(f2 - n1*n2*n3*(n4 - 1)*eps2));
}
