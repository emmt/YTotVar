/*
 * ytotvar.c --
 *
 * Implementation of Yorick interface to relaxed Total Variation (TV).
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>,
 *                          Loïc Denis <loic.denis@univ-st-etienne.fr>,
 *                          Ferréol Soulez <ferreol.soulez@univ-lyon1.fr>
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

#include "totvar.h"

#include <yapi.h>

#ifndef NULL
# define NULL ((void *)0)
#endif

extern void Y_rgl_totvar(int nargs);

void Y_rgl_totvar(int nargs)
{
  double w1, w2, w3;
  double threshold, fx;
  double *x, *gx, *w;
  long dims[Y_DIMSIZE], index, nw, n1, n2, n3, gx_index;
  unsigned int flags;
  int iarg, rank;
  static long options_index = -1L;
  static long threshold_index = -1L;
  static long weight_index = -1L;


  if (options_index < 0L) options_index = yget_global("options", 0);
  if (threshold_index < 0L) threshold_index = yget_global("threshold", 0);
  if (weight_index < 0L) weight_index = yget_global("weight", 0);

  w = NULL;
  x = NULL;
  nw = -1;
  flags = 0;
  rank = 0;
  w1 = w2 = w3 = 1.0;
  threshold = 1E-9;
  n1 = n2 = n3 = 0;
  gx_index = -1L;

  for (iarg = nargs - 1; iarg >= 0; --iarg) {
    index = yarg_key(iarg);
    if (index >= 0) {
      --iarg;
      if (index == options_index) {
        if (! yarg_nil(iarg)) {
          flags = (unsigned int)ygets_l(iarg);
        }
      } else if (index == threshold_index) {
        if (! yarg_nil(iarg)) {
          threshold = ygets_d(iarg);
        }
      } else if (index == weight_index) {
        if (! yarg_nil(iarg)) {
          w = ygeta_d(iarg, NULL, dims);
          if (dims[0] == 0) {
            nw = 0;
          } else if (dims[0] == 1) {
            nw = dims[1];
          }
        }
      } else {
        y_error("unsupported keyword");
      }
    } else if (x == NULL) {
      x = ygeta_d(iarg, NULL, dims);
      rank = dims[0];
      if (rank == 2) {
        n1 = dims[1];
        n2 = dims[2];
      } else if (rank == 3) {
        n1 = dims[1];
        n2 = dims[2];
        n3 = dims[3];
      } else {
        y_error("array X must have 2 or 3 dimensions");
      }
    } else if (gx_index < 0L) {
      gx_index =  yget_ref(iarg);
      if (gx_index < 0L) {
        y_error("expecting simple variable reference for GX");
      }
    } else {
      y_error("too many arguments");
    }
  }
  if (x == NULL){
    y_error("too few arguments");
  }
  if ((rank == 2 && (flags & RGL_TOTVAR_SEPARABLE) == RGL_TOTVAR_SEPARABLE)
      ? (threshold < 0.0) : (threshold <= 0.0)) {
    y_error("bad value for THRESHOLD");
  }
  if (w != NULL) {
    if (nw < 0 || (nw != 0 && nw != rank)) {
      y_error("WEIGHT must be a scalar or a vector of length the number of dimensions of X");
    }
    if (rank == 2) {
      if (nw == 0) {
        w1 = w2 = w[0];
      } else {
        w1 = w[0];
        w2 = w[1];
      }
    } else /* rank = 3 */ {
      if (nw == 0) {
        w1 = w2 = w3 = w[0];
      } else {
        w1 = w[0];
        w2 = w[1];
        w3 = w[2];
      }
    }
    if (w1 < 0.0 || w2 < 0.0 || w3 < 0.0) {
      y_error("weights must all be non-negative");
    }
  }
  if (gx_index >= 0L) {
    /* FIXME: allow to really integrate in pre-existing array GX */
    flags |= RGL_INTEGRATE_GRADIENT; /* this flag is to avoid filling GX with
                                        zeroes twice */
    dims[0] = rank;
    dims[1] = n1;
    dims[2] = n2;
    dims[3] = n3;
    gx = ypush_d(dims);
  } else {
    gx = NULL;
  }
  if (rank == 2) {
    fx = rgl_tv2d(x, n1, n2, w1, w2, threshold, gx, flags);
  } else  /* rank = 3 */ {
    fx = rgl_tv3d(x, n1, n2, n3, w1, w2, w3, threshold, gx, flags);
  }
  if (fx == -1.0) {
    y_error("bug in rgl_totvar");
  }
  if (gx_index >= 0L) {
    yput_global(gx_index, 0);
  }
  ypush_double(fx);
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 78
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
