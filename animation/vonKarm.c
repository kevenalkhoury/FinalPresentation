#include "navier-stokes/centered.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
double dx;

int main() {
//The domain is eight units long, centered vertically.
	L0 = 8.;
	origin (-0.5, -L0/2.);
	N = 512;
	dx=L0/N;
//set a constant viscosity corresponding to a Reynolds number of 160, based on the cylinder diameter (0.125) and the inflow velocity
//set the initial velocity field and tracer concentration
	const face vector muc[] = {0.00078125,0.00078125};
	mu = muc;
	run(); 
}
//fluid is injected on the left boundary with a unit velocity
//tracer is injected in the lower-half of the left boundary. 
//An outflow condition is used on the right boundary
u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
//add a new boundary condition for the cylinder:
//the tangential velocity on the cylinder is set to zero.
bid cylinder;
u.t[cylinder] = dirichlet(0.);

event init (t = 0) {
//make a long channel, 
//top boundary for y>0.5 and the bottom boundary for y<−0.5. 
//The cylinder has a radius of 0.0625
	 mask (y >  0.5 ? top :
	       y < -0.5 ? bottom :
	       sq(x) + sq(y) < sq(0.0625) ? cylinder :
	       none);
//set the initial velocity field
	 foreach()
		      u.x[] = 1.;
}
//check the number of iterations of the Poisson and viscous problems.
event logfile (i++)
	  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
//produce animations of the vorticity and tracer fields
event movies (i += 4; t <= 15.) {
	static FILE * fp = popen ("ppm2gif >vort.ppm", "w");
	scalar vorticity[];
	foreach()
		vorticity[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*dx);
	boundary ({vorticity});
	output_ppm (vorticity, fp, box = {{-0.5,-0.5},{7.5,0.5}},
	min = -10, max = 10, linear = true);
	static FILE * fp1 = popen ("ppm2gif > f.ppm", "w");
	output_ppm (f, fp1, box = {{-0.5,-0.5},{7.5,0.5}},
	linear = true, min = 0, max = 1);
}
//adapt according to the error on the velocity and tracer fields
event adapt (i++) {
	  adapt_wavelet ({u,f}, (double[]){3e-2,3e-2,3e-2}, 9, 4);
}
