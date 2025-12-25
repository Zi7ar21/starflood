#include "solver.h"

#include <math.h>

#include "common.h"
#include "config.h"
#include "types.h"

int solver_run(real* restrict acc, real* restrict pot, real* restrict rho, real* restrict prs, const real* restrict mas, const real* restrict rad, const real* restrict pos, const real* restrict vel, unsigned int N, unsigned int step_number) {
	#ifdef PAIRWISE_SOLVER_DECIMATION
	unsigned int j_length = N / (unsigned int)PAIRWISE_SOLVER_DECIMATION;
	unsigned int j_offset = (step_number % (unsigned int)PAIRWISE_SOLVER_DECIMATION) * j_length;
	#endif

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOADING
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
		real pot_sum = (real)0.0;
		#ifdef ENABLE_KAHAN_SUMMATION
		real acc_c[3] = {(real)0.0, (real)0.0, (real)0.0};
		real pot_c = (real)0.0;
		#endif

		#ifdef ENABLE_SPH
		real rho_sum = (real)0.0;
		#ifdef ENABLE_KAHAN_SUMMATION
		real rho_c = (real)0.0;
		#endif
		#endif

		real m_i = mas[i]; // body i's mass

		// body i's position
		real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		#ifndef PAIRWISE_SOLVER_DECIMATION
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

			// vector from body j's position to body i's position
			real r_ij[3] = {
				r_i[0u] - r_j[0u],
				r_i[1u] - r_j[1u],
				r_i[2u] - r_j[2u]
			};

			// squared distance to body j
			real r2 = (r_ij[0u]*r_ij[0u])+(r_ij[1u]*r_ij[1u])+(r_ij[2u]*r_ij[2u]);

			#ifdef EPSILON
			real inv_r2 = (real)1.0 / ( r2 + (real)(EPSILON*EPSILON) );
				#ifdef STARFLOOD_DOUBLE_PRECISION
				real inv_r1 = (real)1.0 / sqrt( r2 + (real)(EPSILON*EPSILON) );
				#else
				real inv_r1 = (real)1.0 / sqrtf( r2 + (real)(EPSILON*EPSILON) );
				#endif
			#else
			real inv_r2 = (real)0.0 < r2 ? (real)1.0 / r2 : (real)0.0;
				#ifdef STARFLOOD_DOUBLE_PRECISION
				real inv_r1 = (real)0.0 < r1 ? (real)1.0 / sqrt(r1) : (real)0.0;
				#else
				real inv_r1 = (real)0.0 < r1 ? (real)1.0 / sqrtf(r1) : (real)0.0;
				#endif
			#endif

			// gravitational potential of body j
			// (G and m_i are taken into account later)
			real pot_j = m_j * inv_r1;

			#ifdef ENABLE_SPH
			real a = (real)1.0 / (h_j * (real)SQRT_PI);

			#ifdef STARFLOOD_DOUBLE_PRECISION
			real W = (a * a * a) * (real)exp( -r2 / (h_j * h_j) );
			#else
			real W = (a * a * a) * (real)expf( -r2 / (h_j * h_j) );
			#endif

			real rho_j = m_j * W;
			#endif

			real F[3] = {
				pot_j * r_ij[0u] * inv_r2,
				pot_j * r_ij[1u] * inv_r2,
				pot_j * r_ij[2u] * inv_r2
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = F[k] - acc_c[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_c[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += F[k];
				#endif
			}

			{
				#ifdef ENABLE_KAHAN_SUMMATION
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = pot_j - pot_c;
				volatile real t = pot_sum + y;
				volatile real z = t - pot_sum;
				pot_c = z - y;
				pot_sum = t;
				#else
				// Naïve summation
				pot_sum += pot_j;
				#endif
			}

			#ifdef ENABLE_SPH
			{
				#ifdef ENABLE_KAHAN_SUMMATION
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = rho_j - rho_c;
				volatile real t = rho_sum + y;
				volatile real z = t - rho_sum;
				rho_c = z - y;
				rho_sum = t;
				#else
				// Naïve summation
				rho_sum += pot_j;
				#endif
			}
			#endif
		}

		#ifdef PAIRWISE_SOLVER_DECIMATION
		acc[3u*i+0u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[0u];
		acc[3u*i+1u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[1u];
		acc[3u*i+2u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[2u];

		pot[i] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * pot_sum;
		#else
		acc[3u*i+0u] = (real)(-G) * m_i * acc_sum[0u];
		acc[3u*i+1u] = (real)(-G) * m_i * acc_sum[1u];
		acc[3u*i+2u] = (real)(-G) * m_i * acc_sum[2u];

		pot[i] = (real)(-G) * pot_sum;
		#endif

		#ifdef ENABLE_SPH
		rho[i] = rho_sum;
		#endif
	}

	#ifdef ENABLE_SPH
	for(unsigned int i = 0u; i < N; i++) {
		const real k = (real)0.100; // equation of state constant
		const real n = (real)1.000; // polytropic index

		#ifdef STARFLOOD_DOUBLE_PRECISION
		prs[i] = (real)pow(k * rho[i], (real)1.0+(real)1.0/n);
		#else
		prs[i] = (real)powf(k * rho[i], (real)1.0+(real)1.0/n);
		#endif
	}

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOADING
		#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N]) map(to: mas[:N], rad[:N], pos[:3u*N], rho[:N], prs[:N]) 
		#else
		#pragma omp parallel for schedule(dynamic, 128)
		#endif
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real acc_sum[3] = {(real)0.0, (real)0.0, (real)0.0};
		#ifdef ENABLE_KAHAN_SUMMATION
		real acc_c[3] = {(real)0.0, (real)0.0, (real)0.0};
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

		for(unsigned int j = 0u; j < N; j++) {
			if(i == j) {
				continue;
			}

			real m_j = mas[j]; // body j's mass

			real h_j = rad[j]; // body j's smoothing radius

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
				r_i[0u] - r_j[0u],
				r_i[1u] - r_j[1u],
				r_i[2u] - r_j[2u]
			};

			// squared distance to body j
			real r2 = (r_ij[0u]*r_ij[0u])+(r_ij[1u]*r_ij[1u])+(r_ij[2u]*r_ij[2u]);

			if( r2 > (real)(1.000 * 1.000) ) {
				continue;
			}

			// distance to body j
			#ifdef STARFLOOD_DOUBLE_PRECISION
			real r1 = (real)sqrt(r2);
			#else
			real r1 = (real)sqrtf(r2);
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

			#ifdef STARFLOOD_DOUBLE_PRECISION
			real n = (real)(-2.0) * (real)exp( -r2 / (h_j * h_j) ) / (h_j * h_j * h_j * h_j * h_j) / (real)pow( (real)PI, (real)(3.0/2.0) );
			#else
			real n = (real)(-2.0) * expf( -r2 / (h_j * h_j) ) / (h_j * h_j * h_j * h_j * h_j) / (real)powf( (real)PI, (real)(3.0/2.0) );
			#endif

			real W_grad[3] = {
				n * r_ij[0u],
				n * r_ij[1u],
				n * r_ij[2u],
			};

			real F[3] = {
				( prs_i / (rho_i * rho_i + (real)1.0e-6) + prs_j / (rho_j * rho_j + (real)1.0e-6) ) * W_grad[0u],
				( prs_i / (rho_i * rho_i + (real)1.0e-6) + prs_j / (rho_j * rho_j + (real)1.0e-6) ) * W_grad[1u],
				( prs_i / (rho_i * rho_i + (real)1.0e-6) + prs_j / (rho_j * rho_j + (real)1.0e-6) ) * W_grad[2u]
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = F[k] - acc_c[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_c[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += F[k];
				#endif
			}
		}

		acc[3u*i+0u] += -m_i * acc_sum[0u];
		acc[3u*i+1u] += -m_i * acc_sum[1u];
		acc[3u*i+2u] += -m_i * acc_sum[2u];
	}
	#endif

	return STARFLOOD_SUCCESS;
}
