// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "solver.h"

#include <math.h>
#include <stdio.h>

#include "common.h"
#include "config.h"
#include "timing.h"
#include "types.h"

int solver_run(real* restrict acc, real* restrict pot, real* restrict rho, real* restrict prs, const real* restrict mas, const real* restrict rad, const real* restrict pos, const real* restrict vel, unsigned int N, unsigned int step_number) {
	TIMING_INIT();

	#ifdef PAIRWISE_SOLVER_DECIMATION
	unsigned int j_length = N / (unsigned int)PAIRWISE_SOLVER_DECIMATION;
	unsigned int j_offset = (step_number % (unsigned int)PAIRWISE_SOLVER_DECIMATION) * j_length;
	#endif

	TIMING_START();

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_SIM
			#ifdef ENABLE_SPH
			#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N], pot[:N], rho[:N]) map(to: mas[:N], rad[:N], pos[:3u*N]) 
			#else
			#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N], pot[:N]) map(to: mas[:N], pos[:3u*N]) 
			#endif
		#else
		#pragma omp parallel for schedule(dynamic, 128)
		#endif
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real acc_sum[3] = {(real)0.0, (real)0.0, (real)0.0};

		#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
		real acc_com[3] = {(real)0.0, (real)0.0, (real)0.0};
		#endif

		real pot_sum = (real)0.0;

		#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
		real pot_com = (real)0.0;
		#endif

		#ifdef ENABLE_SPH
		real rho_sum = (real)0.0;
		#ifdef ENABLE_KAHAN_SUMMATION
		real rho_com = (real)0.0;
		#endif
		#endif

		real m_i = mas[i]; // body i's mass

		// body i's position
		real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		#ifdef PAIRWISE_SOLVER_DECIMATION
		for(unsigned int j = j_offset; j < (j_offset+j_length); j++) {
		#else
		for(unsigned int j = 0u; j < N; j++) {
		#endif
			/*
			if(i == j) {
				continue;
			}
			*/

			real m_j = mas[j]; // body j's mass

			#ifdef ENABLE_SPH
			real h_j = rad[j]; // body j's smoothing radius
			#endif

			// body j's position
			real r_j[3] = {
				pos[3u*j+0u],
				pos[3u*j+1u],
				pos[3u*j+2u]
			};

			// vector from body j's position to body i's position
			real r_ij[3] = {
				r_i[0] - r_j[0],
				r_i[1] - r_j[1],
				r_i[2] - r_j[2]
			};

			// squared distance to body j
			real r2 = (r_ij[0]*r_ij[0])+(r_ij[1]*r_ij[1])+(r_ij[2]*r_ij[2]);

			#ifdef EPSILON
			real inv_r2 = (real)1.0 /          ( r2 + (real)(EPSILON*EPSILON) );
			real inv_r1 = (real)1.0 / real_sqrt( r2 + (real)(EPSILON*EPSILON) );
			#else
			real inv_r2 = (real)0.0 < r2 ? (real)1.0 /          (r2) : (real)0.0;
			real inv_r1 = (real)0.0 < r2 ? (real)1.0 / real_sqrt(r2) : (real)0.0;
			#endif

			// gravitational potential of body j
			// (G and m_i are taken into account later)
			real pot_j = m_j * inv_r1;

			#ifdef ENABLE_SPH
			real a = (real)1.0 / (h_j * (real)SQRT_PI);

			real W = (a * a * a) * real_exp( -r2 / (h_j * h_j) );

			real rho_j = m_j * W;
			#endif

			real F[3] = {
				pot_j * inv_r2 * r_ij[0],
				pot_j * inv_r2 * r_ij[1],
				pot_j * inv_r2 * r_ij[2]
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = F[k] - acc_com[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_com[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += F[k];
				#endif
			}

			#ifdef ENABLE_SPH
			{
				#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = rho_j - rho_com;
				volatile real t = rho_sum + y;
				volatile real z = t - rho_sum;
				rho_com = z - y;
				rho_sum = t;
				#else
				// Naïve summation
				rho_sum += pot_j;
				#endif
			}
			#endif

			{
				#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = pot_j - pot_com;
				volatile real t = pot_sum + y;
				volatile real z = t - pot_sum;
				pot_com = z - y;
				pot_sum = t;
				#else
				// Naïve summation
				pot_sum += pot_j;
				#endif
			}
		}

		#ifdef PAIRWISE_SOLVER_DECIMATION
		acc[3u*i+0u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[0u];
		acc[3u*i+1u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[1u];
		acc[3u*i+2u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[2u];
		pot[   i   ] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * pot_sum;
		#else
		acc[3u*i+0u] = (real)(-G) * m_i * acc_sum[0u];
		acc[3u*i+1u] = (real)(-G) * m_i * acc_sum[1u];
		acc[3u*i+2u] = (real)(-G) * m_i * acc_sum[2u];
		pot[   i   ] = (real)(-G) * pot_sum;
		#endif

		#ifdef ENABLE_SPH
		rho[   i   ] = rho_sum;
		#endif
	}

	TIMING_STOP();
	TIMING_PRINT("solver_run()", "brute_force");

	#ifdef ENABLE_SPH
	// https://en.wikipedia.org/wiki/Polytrope
	const real eos_k = (real)0.100; // constant of proportionality
	const real eos_n = (real)1.000; // polytropic index

	{
		const real eos_exponent = (real)( (1.0 / (double)eos_n) + 1.0 );

		TIMING_START();

		for(unsigned int i = 0u; i < N; i++) {
			prs[i] = real_pow(eos_k * rho[i], eos_exponent);
		}
	}

	TIMING_STOP();
	TIMING_PRINT("solver_run()", "sph_pressure");
	TIMING_START();

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_SIM
		#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N]) map(to: mas[:N], rad[:N], pos[:3u*N], rho[:N], prs[:N]) 
		#else
		#pragma omp parallel for schedule(dynamic, 128)
		#endif
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real acc_sum[3] = {(real)0.0, (real)0.0, (real)0.0};
		#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
		real acc_com[3] = {(real)0.0, (real)0.0, (real)0.0};
		#endif

		real m_i = mas[i]; // body i's mass

		// body i's position
		real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		real rho_i = rho[i];
		real prs_i = prs[i];

		real rho_i2 = rho_i * rho_i;

		// particle i's contribution to the SPH conservation of momentum equation
		real con_i = (real)0.0 < rho_i2 ? prs_i / rho_i2 : (real)0.0;

		for(unsigned int j = 0u; j < N; j++) {
			/*
			if(i == j) {
				continue;
			}
			*/

			real m_j = mas[j]; // body j's mass

			// body j's position
			real r_j[3] = {
				pos[3u*j+0u],
				pos[3u*j+1u],
				pos[3u*j+2u]
			};

			real rho_j = rho[j];
			real prs_j = prs[j];

			// vector from body j's position to body i's position
			real r_ij[3] = {
				r_i[0] - r_j[0],
				r_i[1] - r_j[1],
				r_i[2] - r_j[2]
			};

			real rho_j2 = rho_j * rho_j;

			// particle j's contribution to the SPH conservation of momentum equation
			real con_j = (real)0.0 < rho_j2 ? prs_j / rho_j2 : (real)0.0;

			// squared distance to body j
			real r2 = (r_ij[0]*r_ij[0])+(r_ij[1]*r_ij[1])+(r_ij[2]*r_ij[2]);

			if( r2 > (real)(1.000 * 1.000) ) {
				continue;
			}

			// distance to body j
			real r1 = real_sqrt(r2);

			#ifdef EPSILON
			real inv_r2 = (real)1.0 /          ( r2 + (real)(EPSILON*EPSILON) );
			real inv_r1 = (real)1.0 / real_sqrt( r2 + (real)(EPSILON*EPSILON) );
			#else
			real inv_r2 = (real)1.0 / r2;
			real inv_r1 = (real)1.0 / r1;
			#endif

			real n = (real)-2.0 * real_exp( -r2 / (h_j * h_j) ) / (INV_POW_PI_3_2 * h_j * h_j * h_j * h_j * h_j);

			real W_grad[3] = {
				n * r_ij[0],
				n * r_ij[1],
				n * r_ij[2],
			};

			real X[3] = {
				(con_i + con_j) * W_grad[0],
				(con_i + con_j) * W_grad[1],
				(con_i + con_j) * W_grad[2]
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = X[k] - acc_com[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_com[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += X[k];
				#endif
			}
		}

		acc[3u*i+0u] += -m_i * acc_sum[0u];
		acc[3u*i+1u] += -m_i * acc_sum[1u];
		acc[3u*i+2u] += -m_i * acc_sum[2u];
	}

	TIMING_STOP();
	TIMING_PRINT("solver_run()", "sph_conservation");
	#endif

	return STARFLOOD_SUCCESS;
}
