/*
 * totvar.i --
 *
 * Yorick interface to relaxed Total Variation (TV).
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

if (is_func(plug_in)) plug_in, "totvar";

local RGL_TOTVAR_ISOTROPIC, RGL_TOTVAR_FORWARD;
local RGL_TOTVAR_ALTERNATIVE, RGL_TOTVAR_SEPARABLE;
extern rgl_totvar;
/* DOCUMENT rgl_totvar(x)
         or rgl_totvar(x, gx)

     This function returns the Total Variation (TV) cost for the, 2-D, 3-D or
     4-D, array X.  Optional argument GX is an output variable used to store
     the partial derivatives of the cost function; if GX is not specified, the
     partial derivatives are not computed.

     The cost is (the sum over every locations of):

        TV = sqrt(NORM2GRAD + THRESHOLD^2) - THRESHOLD

     where NORM2GRAD is (an approximation of) the squared norm of the spatial
     gradient of X along its dimensions.  For instance in 2-D:

        NORM2GRAD(i1,i2)  = WEIGHT(1)*Q1(i1,i2) + WEIGHT(2)*Q2(i1,i2)

     with:

        Q1(i1,i2) = (   ( X(i1+1,i2)   - X(i1,i2)   )^2
                      + ( X(i1+1,i2+1) - X(i1,i2+1) )^2 )/2.0

        Q2(i1,i2) = (   ( X(i1,i2+1)   - X(i1,i2)   )^2
                      + ( X(i1+1,i2+1) - X(i1+1,i2) )^2 )/2.0

     Keyword WEIGHT can be set with the weights along the dimensions of X.
     Set WEIGHT to a scalar to have the same weight along all dimensions.  By
     default, WEIGHT is one for all dimensions.

     Keyword THRESHOLD can be used to set the relaxation parameter.  THRESHOLD
     must be a strictly positive value (to avoid singularities).  This
     parameter can be used to achieve edge-perserving smoothness.

     Keyword OPTIONS can be set with an integer to indicate the TV formula
     to apply:

       RGL_TOTVAR_ISOTROPIC   - use an isotropic definition of TV (this is
                                the default);
       RGL_TOTVAR_FORWARD     - use forward differences for the gradient
                                along the dimensions of X;
       RGL_TOTVAR_ALTERNATIVE - use an alternative definition of isotropic TV
                                which involves the finite differences along
                                the diagonals (only implemented for 2-D arrays);
       RGL_TOTVAR_SEPARABLE   - use a separable definition of TV (only
                                implemented for 2-D arrays);
*/
RGL_TOTVAR_ISOTROPIC    = 0;
RGL_TOTVAR_FORWARD      = 1;
RGL_TOTVAR_ALTERNATIVE  = 2;
RGL_TOTVAR_SEPARABLE    = 3;

extern _rgl_mixed_2dpt;
/* PROTOTYPE
   double mixed_regul_2dpt(double, double, double, double,
                           long, long, long,
                           double array, pointer, int); */

extern _rgl_mixed_3dpt;
/* PROTOTYPE
   double mixed_regul_3dpt(double, double, double, double,
                           long, long, long, long,
                           double array, pointer, int); */

func rgl_mixed_ndpt(args)
/* DOCUMENT rgl_mixed_ndpt(mu1, eps1, mu2, eps2, x);
         or rgl_mixed_ndpt(mu1, eps1, mu2, eps2, x, g, clr);

   SEE ALSO:
 */
{
  local mu1, eps1, mu2, eps2, x, g, clr;
  nargs = args(0);
  if (nargs < 5 || nargs > 7) {
    error, "bad number of arguments";
  }
  eq_nocopy, mu1,  args(1);
  eq_nocopy, eps1, args(2);
  eq_nocopy, mu2,  args(3);
  eq_nocopy, eps2, args(4);
  eq_nocopy, x,    args(5);
  dims = dimsof(x);
  rank = numberof(dims) - 1;
  if (rank < 3 || rank > 4 || structof(x) != double) {
    error, "expecting 3-D or 4-D array of double precision reals";
  }
  if (mu1 < 0 || eps1 <= 0 || mu2 < 0 || eps2 <=0) {
    error, "illegal parameters";
  }
  clr = (nargs < 7 ? 0n : !(!args(7)));
  if (nargs >= 6) {
    eq_nocopy, g, args(6);
    gdims = dimsof(g);
    if (structof(g) != double || numberof(gdims) != numberof(dims)
        || anyof(gdims != dims)) {
      g = array(double, dims);
      args, 6, g;
      clr = 0n;
    }
  }
  if (rank == 3) {
    return _rgl_mixed_2dpt(mu1, eps1, mu2, eps2,
                           dims(2), dims(3), dims(4),
                           x, &g, 0n);
  } else if (rank == 4) {
    return _rgl_mixed_3dpt(mu1, eps1, mu2, eps2,
                           dims(2), dims(3), dims(4), dims(5),
                           x, &g, 0n);
  }
}
wrap_args, rgl_mixed_ndpt;

