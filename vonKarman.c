#include "navier-stokes/centered.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
 const face vector muc[] = {0.00078125,0.00078125};
  mu = muc;
  run(); 
}
u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
bid cylinder;
u.t[cylinder] = dirichlet(0.);

event init (t = 0) {
 mask (y >  0.5 ? top :
	y < -0.5 ? bottom :
	sq(x) + sq(y) < sq(0.0625) ? cylinder :
	none);
 foreach()
    u.x[] = 1.;
}
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
event movies (i += 4; t <= 15.) {
  scalar omega[];
  vorticity (u,omega);
  output_ppm (omega, file = "vort.gif", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true);
  output_ppm (f, file = "f.gif", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = true, min = 0, max = 1);
}
#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D -s ../karman.gfv", "w");
  output_gfs (fp);
}
#endif


event adapt (i++) {
  adapt_wavelet ({u,f}, (double[]){3e-2,3e-2,3e-2}, 9, 4);
}
