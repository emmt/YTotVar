/*
 * totvar.h --
 *
 * Definitions for relaxed Total Variation (TV).
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

#ifndef _TOTVAR_H
#define _TOTVAR_H

/* Options. */
#define RGL_TOTVAR_ISOTROPIC    (0)
#define RGL_TOTVAR_FORWARD      (1)
#define RGL_TOTVAR_ALTERNATIVE  (2)
#define RGL_TOTVAR_SEPARABLE    (3)

/* Flags to combine with options. */
#define RGL_NO_GRADIENT         (0)       /* do not compute gradient */
#define RGL_STORE_GRADIENT      (1 << 3)  /* compute gradient, pre-fill
                                             gradient array with zeros */
#define RGL_INTEGRATE_GRADIENT  (1 << 4)  /* integrate gradient */


#ifdef __cplusplus
extern "C" {
#endif

extern double rgl_tv2d(const double x[],
                       const long n1, const long n2,
                       const double w1, const double w2,
                       const double eps, double gx[],
                       unsigned int flags);

extern double rgl_tv3d(const double x[],
                       const long n1, const long n2, const long n3,
                       const double w1, const double w2, const double w3,
                       const double eps, double gx[],
                       unsigned int flags);

extern double rgl_tv4d(const double x[],
                       const long n1, const long n2,
                       const long n3, const long n4,
                       const double w1, const double w2,
                       const double w3, const double w4,
                       const double eps, double gx[],
                       unsigned int flags);

#ifdef __cplusplus
}
#endif

#endif /* _TOTVAR_H */
