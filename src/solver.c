#include "solver.h"

#include <math.h>

#include "common.h"
#include "config.h"
#include "types.h"

int solver_run(real* volatile pot, real* volatile acc, const real* volatile mas, const real* volatile pos, unsigned int N, unsigned int step_number) {
	#ifdef N_DIV
	unsigned int j_length = N / (unsigned int)N_DIV;
	unsigned int j_offset = (step_number % (unsigned int)N_DIV) * j_length;
	#endif

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOADING
		#pragma omp target teams distribute parallel for map(to: mas[:N], pos[:3u*N]) map(tofrom: pot[:N], acc[:3u*N]) 
		#else
		#pragma omp parallel for schedule(dynamic, 128)
		#endif
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real U_sum = (real)0.0;

		#ifdef SOLVER_USE_KAHAN_SUMMATION_ENERGY
		real U_c = (real)0.0;
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

		real m_i = mas[i]; // Body mass

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

			real m_j = mas[j];

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

			#ifdef EPSILON
			real inv_r2 = (real)1.0 / (       r2+(real)(EPSILON*EPSILON)  );
			real inv_r1 = (real)1.0 / ( sqrtf(r2+(real)(EPSILON*EPSILON)) );
			#else
			real inv_r2 = (real)1.0 / (       r2  );
			real inv_r1 = (real)1.0 / ( sqrtf(r2) );
			#endif

			// Gravitational potential
			real U_ij = -(real)G * m_i * m_j * inv_r2;

			{
				#ifdef SOLVER_USE_KAHAN_SUMMATION_ENERGY
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = U_ij - U_c;
				volatile real t = U_sum + y;
				volatile real z = t - U_sum;
				U_c = z - y;
				U_sum = t;
				#else
				// Naïve summation
				U_sum += U_ij;
				#endif
			}

			// Force acting upon body i by body j
			real F_ij[3] = {
				U_ij * r_ij[0u] * inv_r1,
				U_ij * r_ij[1u] * inv_r1,
				U_ij * r_ij[2u] * inv_r1
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
		}

		pot[i] = U_sum;

		#ifndef N_DIV
		acc[3u*i+0u] = F_sum[0u];
		acc[3u*i+1u] = F_sum[1u];
		acc[3u*i+2u] = F_sum[2u];
		#else
		acc[3u*i+0u] = (real)N_DIV * F_sum[0u];
		acc[3u*i+1u] = (real)N_DIV * F_sum[1u];
		acc[3u*i+2u] = (real)N_DIV * F_sum[2u];
		#endif
	}

	return STARFLOOD_SUCCESS;
}
