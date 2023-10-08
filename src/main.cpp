#include <cstdio>
#include <cstdlib>
#include <cstdint>

//#include "rng.h"

// OpenMP
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define     PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define INV_PI 0.3183098861837906715377675267450287240689192914809128974953346881

#define     TAU 6.2831853071795864769252867665590057683943387987502116419498891846
#define INV_TAU 0.1591549430918953357688837633725143620344596457404564487476673441

double gaussian2d(double x0, double x1, double sigma) {
	sigma = sigma > 0.0 ? sigma : 1.0;
	double sigma0 = 1.0 / sigma;
	double sigma1 = 1.0 / (2.0*sigma*sigma);
	return sigma0*INV_TAU*exp(-sigma1*((x0*x0)+(x1*x1)));
}

// Triple32 Hash: https://nullprogram.com/blog/2018/07/31/
// this hash in particular is remarkable because it's a statistically perfect 32-bit integer hash (of this kind)
inline uint32_t triple32(uint32_t x) {
	x ^= x >> (uint32_t)17u;
	x *= (uint32_t)0xED5AD4BBu;
	x ^= x >> (uint32_t)11u;
	x *= (uint32_t)0xAC4C1B51u;
	x ^= x >> (uint32_t)15u;
	x *= (uint32_t)0x31848BABu;
	x ^= x >> (uint32_t)14u;
	return x;
}

float urand(uint32_t* ns) {
	uint32_t _ns = *ns;

	_ns = triple32(_ns);

	*ns = _ns;

	return (2.32830629776081821092e-10f)*(float)_ns;
}
#define N 1000

typedef double real;

real W(real x, real y, real z, real h) {
	/*
	Gausssian Smoothing kernel (3D)
	x     is a vector/matrix of x positions
	y     is a vector/matrix of y positions
	z     is a vector/matrix of z positions
	h     is the smoothing length
	w     is the evaluated smoothing function
	*/

	real r = sqrt(x*x+y*y+z*z);

	real w = pow((1.0/(h*sqrt(PI))),3.0) * exp( -(r*r)/ (h*h));

	return w;
}

void gradW(real x, real y, real z, real h, real* grad) {
	/*
	Gradient of the Gausssian Smoothing kernel (3D)
	x     is a vector/matrix of x positions
	y     is a vector/matrix of y positions
	z     is a vector/matrix of z positions
	h     is the smoothing length
	wx, wy, wz     is the evaluated gradient
	*/

	real r = sqrt(x*x+y*y+z*z);

	const real a = -0.0718348488500666246756327934510771021128743773739306898485328149; // -2/(5pi^(3/2))

	real poo = a*exp(-(r*r) / (h*h))*(1.0 / h);
	grad[0] = poo * x;
	grad[1] = poo * y;
	grad[2] = poo * z;

	//return wx, wy, wz;
}

/*
def getPairwiseSeparations( ri, rj ) {
	"""
	Get pairwise desprations between 2 sets of coordinates
	ri    is an M x 3 matrix of positions
	rj    is an N x 3 matrix of positions
	dx, dy, dz   are M x N matrices of separations
	"""
	
	M = ri.shape[0]
	N = rj.shape[0]
	
	# positions ri = (x,y,z)
	rix = ri[:,0].reshape((M,1))
	riy = ri[:,1].reshape((M,1))
	riz = ri[:,2].reshape((M,1))
	
	# other set of points positions rj = (x,y,z)
	rjx = rj[:,0].reshape((N,1))
	rjy = rj[:,1].reshape((N,1))
	rjz = rj[:,2].reshape((N,1))
	
	# matrices that store all pairwise particle separations: r_i - r_j
	dx = rix - rjx.T
	dy = riy - rjy.T
	dz = riz - rjz.T
	
	return dx, dy, dz
}
*/

real getDensity(real* x, real* r, real m, real h) {
	/*
	Get Density at sampling loctions from SPH particle distribution
	r     is an M x 3 matrix of sampling locations
	pos   is an N x 3 matrix of SPH particle positions
	m     is the particle mass
	h     is the smoothing length
	rho   is M x 1 vector of densities
	*/

	/*
	//M = r.shape[0];

	dx, dy, dz = getPairwiseSeparations( r, pos );

	rho = np.sum( m * W(dx, dy, dz, h), 1 ).reshape((M,1));

	real rho = 0;

	for(int i = 0; i < N; i++) {
		rho += m * W(dx, dy, dz, h);
	}

	return rho;
	*/

	double rho = (double)0;

	for(int i = 0; i < N; i++) {
		rho += m * W(r[3*i+0]-x[0], r[3*i+1]-x[1], r[3*i+2]-x[2], h);
	}

	return rho;
}

real getPressure(real rho, real k, real n) {
	/*
	Equation of State
	rho   vector of densities
	k     equation of state constant
	n     polytropic index
	P     pressure
	*/

	real P = k * pow(rho, 1.0+1.0/n);

	return P;
}

void getAcc(real* r, real* v, real* rho, real* P, real m, real h, real k, real n, real lambda, real nu, real* a) {
	/*
	Calculate the acceleration on each SPH particle
	pos   is an N x 3 matrix of positions
	vel   is an N x 3 matrix of velocities
	m     is the particle mass
	h     is the smoothing length
	k     equation of state constant
	n     polytropic index
	lmbda external force constant
	nu    viscosity
	a     is N x 3 matrix of accelerations
	*/

	//N = pos.shape[0];

	// Calculate densities at the position of the particles
	for(int i = 0; i < N; i++) {
		real cr[3] = {r[3*i+0], r[3*i+1], r[3*i+2]};

		rho[i] = getDensity(cr, r, m, h);
	}

	// Get the pressures
	for(int i = 0; i < N; i++) P[i] = getPressure(rho[i], k, n);

	for(int i = 0; i < N; i++) {
		a[3*i+0] = (real)0;
		a[3*i+0] = (real)0;
		a[3*i+0] = (real)0;

		for(int j = 0; j < N; j++) {
			if(i == j) continue;

			real dx = r[3*j+0]-r[3*i+0];
			real dy = r[3*j+1]-r[3*i+1];
			real dz = r[3*j+2]-r[3*i+2];
			real grad[3];
			gradW(dx,dy,dz,h,grad);
			a[3*i+0] += m*((P[i]/(rho[i]*rho[i]))+(P[j]/(rho[j]*rho[j])))*grad[0];
			a[3*i+1] += m*((P[i]/(rho[i]*rho[i]))+(P[j]/(rho[j]*rho[j])))*grad[1];
			a[3*i+2] += m*((P[i]/(rho[i]*rho[i]))+(P[j]/(rho[j]*rho[j])))*grad[2];

			// Get pairwise distances and gradients


			// Add Pressure contribution to accelerations
			//ax = - np.sum( m * ( P/(rho**2) + P.T/rho.T**2  ) * dWx, 1).reshape((N,1));
			//ay = - np.sum( m * ( P/(rho**2) + P.T/rho.T**2  ) * dWy, 1).reshape((N,1));
			//az = - np.sum( m * ( P/(rho**2) + P.T/rho.T**2  ) * dWz, 1).reshape((N,1));

			// pack together the acceleration components
			//a = np.hstack((ax,ay,az));
		}
	}

	// Add external potential force
	for(int i = 0; i < (N * 3); i++) a[i] -= lambda * r[i];

	// Add viscosity
	for(int i = 0; i < (N * 3); i++) a[i] -= nu * v[i];

	//return a;
}

int main(int argc, char** argv) {
	if(argc > 1) {
		printf("%s: Invalid Parameter '%s'\n", argv[0], argv[1]);

		return EXIT_FAILURE;
	}

	//const int N = 400; // number of particles
	const real dt = 0.04; // timestep
	const real M = 2.0; // star mass
	const real R = 0.75; // star radius
	const real h = 0.1; // smoothing length
	const real k = 0.1; // equation of state constant
	const real n = 1.0; // polytropic index
	const real nu = 0.5; // damping

	size_t r_size = (N * 3) * sizeof(real);
	size_t v_size = (N * 3) * sizeof(real);
	size_t a_size = (N * 3) * sizeof(real);
	size_t t_size = (N * 3) * sizeof(real);
	size_t P_size = (N * 1) * sizeof(real);
	size_t rho_size = N * sizeof(real);

	real* r = (real*)malloc(r_size); // position buffer
	real* v = (real*)malloc(v_size); // velocity buffer
	real* a = (real*)malloc(a_size); // acceleration buffer
	real* P = (real*)malloc(P_size); // pressure
	real* rho = (real*)malloc(rho_size); // interpolated density

	for(int i = 0; i < n * 3; i++) a[i] = (real)0;
	// Generate Initial Conditions
	uint32_t ns = (uint32_t)1u;

	real m = (real)((double)M/(double)N); // single particle mass

	m = (real)0.5;

	//real lambda = 2.0*k*(1.0+n)*pow(PI, -3.0/(2.0*n)) * pow((M*gamma(5.0/2.0+n)/(R*R*R)/gamma(1.0+n)),(1.0/n)) / (R*R); // ~ 2.01
	real lambda = 2.0;

	// randomly selected positions and velocities
	for(int i = 0; i < N; i++) {
		ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

		float z0 = urand(&ns);
		float z1 = urand(&ns);
		float z2 = urand(&ns);
		float z3 = urand(&ns);

		r[3*i+0] = 0.1*sqrt(-2.0*log(z0))*cos(TAU*z2);
		r[3*i+1] = 0.1*sqrt(-2.0*log(z0))*sin(TAU*z2);
		r[3*i+2] = 0.1*sqrt(-2.0*log(z1))*cos(TAU*z3);
	}

	for(int i = 0; i < (N * 3); i++) v[i] = (real)0;
	
	// calculate initial gravitational accelerations
	getAcc(r, v, rho, P, m, h, k, n, lambda, nu, a);
	int image_w = 640, image_h = 480;

	double t0, t1;

	t0 = omp_get_wtime();
	t1 = omp_get_wtime();

	size_t image_size = (image_w * image_h) * 3 * sizeof(float);

	float* image = (float*)malloc(image_size);

	stbi_flip_vertically_on_write(1);

	const int num_steps = 600;

	//const real dt = 0.1;

	for(int step_num = 0; step_num < num_steps; step_num++) {
		// render and write
		{
			// clear the image buffer
			//#pragma omp simd
			for(int i = 0; i < (image_w * image_h) * 3; i++) image[i] = (float)0;

			for(int i = 0; i < N; i++) {
				// vec2 uv = (fragCoord - 0.5 * resolution) / resolution.y
				real t = 0.04*(real)step_num;
				real uv[2] = {r[2 * i + 0]*cos(t)+r[2 * i + 2]*sin(t), r[2 * i + 1]};

				uv[0] *= (real)2.0; // 0.05
				uv[1] *= (real)2.0;

				int coord[2] = {
				(int)((real)image_h * uv[0] + (real)0.5 * (real)image_w),
				(int)((real)image_h * uv[1] + (real)0.5 * (real)image_h)};

				if(0 <= coord[0] && coord[0] < image_w && 0 <= coord[1] && coord[1] < image_h) {
					image[3 * (image_w * coord[1] + coord[0]) + 0] += 0.5f;
					image[3 * (image_w * coord[1] + coord[0]) + 1] += 0.5f;
					image[3 * (image_w * coord[1] + coord[0]) + 2] += 0.5f;
				}
			}

			char file_name[FILENAME_MAX];

			snprintf(file_name, FILENAME_MAX, "out/step_%04d.hdr", step_num);

			t0 = omp_get_wtime();

			if(stbi_write_hdr(file_name, image_w, image_h, 3, image) == 0) {
				printf("Failure on step %d\n", step_num);
			}

			t1 = omp_get_wtime();

			//fprintf(logfile, "%.6f,", 1000.0*(t1-t0));
		}

		// kick
		{
			for(int i = 0; i < (N * 3); i++) v[i] += (real)0.5 * dt * a[i];
		}

		// drift
		{
			for(int i = 0; i < (N * 3); i++) r[i] += dt * v[i];
		}

		/*
		update accelerations
				acc = getAcc( pos, vel, m, h, k, n, lmbda, nu )
		*/
		for(int i = 0; i < (N * 3); i++) a[i] = (real)0;
		getAcc(r, v, rho, P, m, h, k, n, lambda, nu, a);
		// update acceleration
		//{
		//	for(int i = 0; i < n * 3; i++) a[i] = (real)0;
		//}

		// kick
		{
			for(int i = 0; i < (N * 3); i++) v[i] += (real)0.5 * dt * a[i];
		}
	}

	free(r);
	free(v);
	free(a);
	free(P);
	free(rho);
	free(image);

	return EXIT_SUCCESS;
}
