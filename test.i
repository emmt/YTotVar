//plug_dir,".";
//include, "./totvar.i";
func rgl_totvar_test(dims, threshold=, weight=, options=)
{
  local gx;
  epsilon = 2.3E-16; /* machine precision */
  x = random(dims) - 0.5;
  f0 = rgl_totvar(x, threshold=threshold, options=options, weight=weight);
  f1 = rgl_totvar(x, gx, threshold=threshold, options=options, weight=weight);
  if (abs(f1 - f0) > epsilon*abs(f0)) {
    write, format="*** ERRORR: F1 - F0 = %+#E  (F0 = %+#E)\n", f1-f0, f0;
  }

  errmax = 0.0;
  tiny = 1E-5;
  xnrm = sqrt(sum(x*x));
  level = tiny*xnrm/numberof(x);
  for (j = 1; j <= 40; ++j) {
    s = level*(random(dimsof(x)) - 0.5);
    f1 = rgl_totvar(x + s, threshold=threshold, options=options, weight=weight);
    sg = sum(s*gx);
    errmax = max(errmax, abs(sg - (f1 - f0)));
    //write, format="s'.g = %+#E    F1 - F0 = %+#E    RATIO = %+#.2E\n",
    //  sg, f1 - f0, (f1 - f0)/sg;
  }
  write, format="GRAD. ERR. MAX. = %+#E\n", errmax/abs(f0)/level;

}

func rgl_runtests
{
  write, format="\n%s\n", "TESTS IN 2D";
  rgl_totvar_test, [2,43,27], threshold = 0.1;
  rgl_totvar_test, [2,43,27], threshold = 0.1, options=RGL_TOTVAR_ISOTROPIC;
  rgl_totvar_test, [2,43,27], threshold = 0.1, options=RGL_TOTVAR_FORWARD;
  rgl_totvar_test, [2,43,27], threshold = 0.1, options=RGL_TOTVAR_ALTERNATIVE;
  rgl_totvar_test, [2,43,27], threshold = 0.1, options=RGL_TOTVAR_SEPARABLE;
  rgl_totvar_test, [2,43,27], threshold = 0.0, options=RGL_TOTVAR_SEPARABLE;

  rgl_totvar_test, [2,43,27], threshold = 0.1, weight=[2,3];
  rgl_totvar_test, [2,43,27], threshold = 0.1, weight=[2,3], options=RGL_TOTVAR_ISOTROPIC;
  rgl_totvar_test, [2,43,27], threshold = 0.1, weight=[2,3], options=RGL_TOTVAR_FORWARD;
  rgl_totvar_test, [2,43,27], threshold = 0.1, weight=[2,3], options=RGL_TOTVAR_ALTERNATIVE;
  rgl_totvar_test, [2,43,27], threshold = 0.1, weight=[2,3], options=RGL_TOTVAR_SEPARABLE;

  write, format="\n%s\n", "TESTS IN 3D";
  rgl_totvar_test, [3,43,27,16], threshold = 0.1;
  rgl_totvar_test, [3,43,27,16], threshold = 0.1, options=RGL_TOTVAR_ISOTROPIC;
  rgl_totvar_test, [3,43,27,16], threshold = 0.1, options=RGL_TOTVAR_FORWARD;
  rgl_totvar_test, [3,43,27,16], threshold = 0.1, weight=[2,3,4];
  rgl_totvar_test, [3,43,27,16], threshold = 0.1, weight=[2,3,4], options=RGL_TOTVAR_ISOTROPIC;
  rgl_totvar_test, [3,43,27,16], threshold = 0.1, weight=[2,3,4], options=RGL_TOTVAR_FORWARD;
}


func __totvar_denoise_f(x, &gx, data)
{
  //local data;
  //extern __totvar_denoise_data;
  //eq_nocopy, data, __totvar_denoise_data;

  gx = 1;
  fx = rgl_totvar(x, gx, threshold=data.threshold,
                  options=data.options, weight=data.weight);
  r = data.y - x;
  fx += data.alpha*sum(r*r);
  gx -= 2.0*data.alpha*r;
  return fx;
}

func totvar_denoise(y, x, mu=, threshold=, weight=, options=,
                    verb=, maxiter=, maxeval=, frtol=, fatol=)
{
  if (! is_func(op_mnb)) include, "OptimPack1.i", 1;
  if (is_void(mu)) mu = 1.0;
  data = h_new(y=double(y), threshold=threshold,
               weight=weight, options=options,
               alpha=1.0/mu);
  if (is_void(x)) x = double(y);
  return op_mnb(__totvar_denoise_f, x, extra=data,
                verb=verb, maxiter=maxiter, maxeval=maxeval,
                frtol=frtol, fatol=fatol);
}

func test_rgl_square(x) { return x*x; }

func test_rgl_mixed(mu1, eps1, mu2, eps2, x, g, clr)
{
  dims = dimsof(x);
  rank = numberof(dims) - 1;
  if (clr && ! is_void(g)) {
    g(*) = 0;
  }
  q = test_rgl_square;
  i = 1:-1;
  j = 2:0;

  /* along leading dimensions */
  f1 = 0.0;
  if (rank == 3) {
    s = sqrt(2);
    mu1  /= s;
    eps1 *= s;
    r = sqrt(q(x(i,i,) - x(i,j,)) +
             q(x(j,i,) - x(j,j,)) +
             q(x(i,i,) - x(j,i,)) +
             q(x(i,j,) - x(j,j,)) + eps1*eps1);
    f1 = mu1*(sum(r) - numberof(r)*eps1);
    if (! is_void(g)) {
      r = mu1/r;
      g(i,i,) += (2*x(i,i,) - x(i,j,) - x(j,i,))*r;
      g(i,j,) += (2*x(i,j,) - x(i,i,) - x(j,j,))*r;
      g(j,i,) += (2*x(j,i,) - x(i,i,) - x(j,j,))*r;
      g(j,j,) += (2*x(j,j,) - x(i,j,) - x(j,i,))*r;
    }
  } else if (rank == 4) {
    s = 2.0;
    mu1  /= s;
    eps1 *= s;
    r = sqrt(q(x(i,i,i,) - x(i,i,j,)) +
             q(x(i,j,i,) - x(i,j,j,)) +
             q(x(j,i,i,) - x(j,i,j,)) +
             q(x(j,j,i,) - x(j,j,j,)) +
             q(x(i,i,i,) - x(i,j,i,)) +
             q(x(i,i,j,) - x(i,j,j,)) +
             q(x(j,i,i,) - x(j,j,i,)) +
             q(x(j,i,j,) - x(j,j,j,)) +
             q(x(i,i,i,) - x(j,i,i,)) +
             q(x(i,i,j,) - x(j,i,j,)) +
             q(x(i,j,i,) - x(j,j,i,)) +
             q(x(i,j,j,) - x(j,j,j,)) + eps1*eps1);
    f1 = mu1*(sum(r) - numberof(r)*eps1);
    if (! is_void(g)) {
      r = mu1/r;
      g(i,i,i,) += (3*x(i,i,i,) - x(i,i,j,) - x(i,j,i,) - x(j,i,i,))*r;
      g(i,i,j,) += (3*x(i,i,j,) - x(i,i,i,) - x(i,j,j,) - x(j,i,j,))*r;
      g(i,j,i,) += (3*x(i,j,i,) - x(i,j,j,) - x(i,i,i,) - x(j,j,i,))*r;
      g(i,j,j,) += (3*x(i,j,j,) - x(i,j,i,) - x(i,i,j,) - x(j,j,j,))*r;
      g(j,i,i,) += (3*x(j,i,i,) - x(j,i,j,) - x(j,j,i,) - x(i,i,i,))*r;
      g(j,i,j,) += (3*x(j,i,j,) - x(j,i,i,) - x(j,j,j,) - x(i,i,j,))*r;
      g(j,j,i,) += (3*x(j,j,i,) - x(j,j,j,) - x(j,i,i,) - x(i,j,i,))*r;
      g(j,j,j,) += (3*x(j,j,j,) - x(j,j,i,) - x(j,i,j,) - x(i,j,j,))*r;
    }
  }

  /* along last dimension */
  t = x(..,i) - x(..,j);
  r = sqrt(t*t + eps2*eps2);
  f2 = mu2*(sum(r) - numberof(r)*eps2);
  if (! is_void(g)) {
    t *= (mu2/r);
    g(..,i) += t;
    g(..,j) -= t;
  }
  return (f1 + f2);
}

func rgl_mixed_runtests(dims, mu1=, eps1=, mu2=, eps2=)
{
  if (is_void(dims)) {
    dims = [3, 7,8,9];
    rgl_mixed_runtests, dims, mu1=mu1, eps1=eps1, mu2=mu2, eps2=eps2;
    dims = [4, 6,7,8,9];
    rgl_mixed_runtests, dims, mu1=mu1, eps1=eps1, mu2=mu2, eps2=eps2;
    return;
  }
  a = random(dims);
  b = random_n(dims);
  if (is_void( mu1))  mu1 = 1.1;
  if (is_void(eps1)) eps1 = 0.3;
  if (is_void( mu2))  mu2 = 2.4;
  if (is_void(eps2)) eps2 = 0.7;

  g1 = array(double, dims);
  g2 = array(double, dims);

  f1 = rgl_mixed_ndpt(mu1, eps1, mu2, eps2, a);
  f2 = test_rgl_mixed(mu1, eps1, mu2, eps2, a);
  write, "f(x) without gradient: ", f1, f2, abs(f1 - f2);

  f1 = rgl_mixed_ndpt(mu1, eps1, mu2, eps2, a, g1,1);
  f2 = test_rgl_mixed(mu1, eps1, mu2, eps2, a, g2,1);
  write, "f(x) with gradient:    ", f1, f2, abs(f1 - f2);
  stat, g1, g2, g1-g2;

  t = 1e-5;
  f3 = rgl_mixed_ndpt(mu1, eps1, mu2, eps2, a + t*b);
  write, "   compiled:", (f3 - f1)/t, sum(b*g1);
  f4 = test_rgl_mixed(mu1, eps1, mu2, eps2, a + t*b);
  write, "interpreted:", (f4 - f2)/t, sum(b*g2);
}
