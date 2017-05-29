 	
%  function m_regul_totvar_2d_complex
% 
%  [fx, gxr, gxi] = m_regul_totvar_2d_complex(xr, xi, w1, w2, eps, flags);
%                                       
%  MexFile to call from Matlab the rgl_tv2d_complex function of the totvar library 
%  which implements Definitions for relaxed Total Variation (TV)
%  for images with complex values.
%  
%  Compilation :
%        mex m_regul_totvar_2d_complex.c totvar.o
%   
%  Mai 2017
%  LaHC
%  Fabien Momey
%   
%  *-----------------------------------------------------------------------------
%  * Information about totvar library
%  *-----------------------------------------------------------------------------
%  *
%  * Copyright (C) 2010-2014, Éric  <eric.thiebaut@univ-lyon1.fr>,
%  *                          Loïc Denis <loic.denis@univ-st-etienne.fr>,
%  *                          Ferréol Soulez <ferreol.soulez@univ-lyon1.fr>
%  *
%  * This software is governed by the CeCILL-C license under French law and
%  * abiding by the rules of distribution of free software.  You can use, modify
%  * and/or redistribute the software under the terms of the CeCILL-C license as
%  * circulated by CEA, CNRS and INRIA at the following URL
%  * "http://www.cecill.info".
%  *
%  * As a counterpart to the access to the source code and rights to copy,
%  * modify and redistribute granted by the license, users are provided only
%  * with a limited warranty and the software's author, the holder of the
%  * economic rights, and the successive licensors have only limited liability.
%  *
%  * In this respect, the user's attention is drawn to the risks associated with
%  * loading, using, modifying and/or developing or reproducing the software by
%  * the user in light of its specific status of free software, that may mean
%  * that it is complicated to manipulate, and that also therefore means that it
%  * is reserved for developers and experienced professionals having in-depth
%  * computer knowledge. Users are therefore encouraged to load and test the
%  * software's suitability as regards their requirements in conditions enabling
%  * the security of their systems and/or data to be ensured and, more
%  * generally, to use and operate it in the same conditions as regards
%  * security.
%  *
%  * The fact that you are presently reading this means that you have had
%  * knowledge of the CeCILL-C license and that you accept its terms.
%  *
%  *-----------------------------------------------------------------------------
