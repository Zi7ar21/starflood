#include "solver.h"

#include <math.h>

#include "common.h"
#include "config.h"
#include "types.h"

int solver_run(real* volatile acc, real* volatile pot, const real* volatile mas, const real* volatile rad, const real* volatile pos, unsigned int N, unsigned int step_number) {
	#ifdef N_DIV
	unsigned int j_length = N / (unsigned int)N_DIV;
	unsigned int j_offset = (step_number % (unsigned int)N_DIV) * j_length;
	#endif

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOADING
			#ifdef ENABLE_SPH
			#pragma omp target teams distribute parallel for map(to: mas[:N], rad[:N], pos[:3u*N]) map(tofrom: acc[:3u*N], pot[:N])
			#else
			#pragma omp target teams distribute parallel for map(to: mas[:N], pos[:3u*N]) map(tofrom: acc[:3u*N], pot[:N])
			#endif
		#else
		#pragma omp parallel for schedule(dynamic, 128)
		#endif
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real V_sum = (real)0.0;

		#ifdef SOLVER_USE_KAHAN_SUMMATION_ENERGY
		real V_c = (real)0.0;
		#endif

		real F_sum[3] = {
			(real)0.0,
			(real)0.0,
			(real)0.0
		};

		#ifdef SOLVER_USE_KAHAN_SUMMATION
		real F_c[3] = {
			(real)0.0,
			(real)0.0,
			(real)0.0
		};
		#endif

		real m_i = mas[i]; // body i's mass

		#ifdef ENABLE_SPH
		real h_i = rad[i]; // body i's smoothing radius
		#endif

		// body i's position
		real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		#ifndef N_DIV
		for(unsigned int j = 0u; j < N; j++) {
		#else
		for(unsigned int j = j_offset; j < (j_offset+j_length); j++) {
		#endif
			if(i == j) {
				continue;
			}

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

			// Vector to position i from position j
			real r_ij[3] = {
				r_i[0u] - r_j[0u],
				r_i[1u] - r_j[1u],
				r_i[2u] - r_j[2u]
			};

			real r2 = (r_ij[0u]*r_ij[0u])+(r_ij[1u]*r_ij[1u])+(r_ij[2u]*r_ij[2u]);

			#ifdef STARFLOOD_DOUBLE_PRECISION
			real r1 = (real)sqrt(r1);
			#else
			real r1 = (real)sqrtf(r1);
			#endif

			#ifdef EPSILON
			real inv_r2 = (real)1.0 / ( r2 + (real)(EPSILON*EPSILON) );
				#ifdef STARFLOOD_DOUBLE_PRECISION
				real inv_r1 = (real)1.0 / sqrt( r2 + (real)(EPSILON*EPSILON) );
				#else
				real inv_r1 = (real)1.0 / sqrtf( r2 + (real)(EPSILON*EPSILON) );
				#endif
			#else
			real inv_r2 = (real)1.0 / r2;
			real inv_r1 = (real)1.0 / r1;
			#endif

			// Gravitational potential
			real U_j = -(real)G * m_j * inv_r1;

			#ifdef ENABLE_SPH
			const real k = (real)1.000; // equation of state constant
			const real n = (real)1.000; // polytropic index

			real P = k * rho (1 + 1/n); // pressure
			real a = 1.0 / (h_i * SQRT_PI);

			real W = (a * a * a) * expf(-r2 / h * h);

			real W_grad[3] = {
			};
			#endif

			// Force acting upon body i by body j
			real F_ij[3] = {
				U_j * r_ij[0u] * inv_r2,
				U_j * r_ij[1u] * inv_r2,
				U_j * r_ij[2u] * inv_r2
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef SOLVER_USE_KAHAN_SUMMATION
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = F_ij[k] - F_c[k];
				volatile real t = F_sum[k] + y;
				volatile real z = t - F_sum[k];
				F_c[k] = z - y;
				F_sum[k] = t;
				#else
				// Naïve summation
				F_sum[k] += F_ij[k];
				#endif
			}

			{
				#ifdef SOLVER_USE_KAHAN_SUMMATION_ENERGY
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = U_j - V_c;
				volatile real t = V_sum + y;
				volatile real z = t - V_sum;
				V_c = z - y;
				V_sum = t;
				#else
				// Naïve summation
				V_sum += V_ij;
				#endif
			}
		}

		#ifndef N_DIV
		acc[3u*i+0u] = m_i * F_sum[0u];
		acc[3u*i+1u] = m_i * F_sum[1u];
		acc[3u*i+2u] = m_i * F_sum[2u];
		#else
		acc[3u*i+0u] = (real)N_DIV * m_i * F_sum[0u];
		acc[3u*i+1u] = (real)N_DIV * m_i * F_sum[1u];
		acc[3u*i+2u] = (real)N_DIV * m_i * F_sum[2u];
		#endif

		pot[i] = m_i * V_sum;
	}

	return STARFLOOD_SUCCESS;
}
