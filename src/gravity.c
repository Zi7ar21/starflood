// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "gravity.h"

#include <math.h>
#include <stdio.h>

#include "common.h"
#include "config.h"
#include "rng.h"
#include "timing.h"
#include "types.h"

int solve_gravity(sim_t* restrict sim_ptr) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	unsigned int N = sim.N;
	real* pos = sim_find(&sim, SIM_POS);
	real* vel = sim_find(&sim, SIM_VEL);
	real* acc = sim_find(&sim, SIM_ACC);
	real* mas = sim_find(&sim, SIM_MAS);
	real* pot = sim_find(&sim, SIM_POT);

	#ifdef ENABLE_GRID
	grid_t grid = sim.grid;

	unsigned int grid_dim[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

	real bounds_min[3] = {grid.bounds_min[0], grid.bounds_min[1], grid.bounds_min[2]};
	real bounds_max[3] = {grid.bounds_max[0], grid.bounds_max[1], grid.bounds_max[2]};

	real bounds_wid[3] = {
		bounds_max[0] - bounds_min[0],
		bounds_max[1] - bounds_min[1],
		bounds_max[2] - bounds_min[2]
	};

	f32* grid_poten = (f32*)grid.mem;
	i32* grid_accum = (i32*)&grid_poten[grid_dim[2] * grid_dim[1] * grid_dim[0]];

	TIMING_START();

	for(unsigned int i = 0u; i < grid.dim[2] * grid.dim[1] * grid.dim[0]; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		grid_accum[i] = (i32)0;
	}

	TIMING_STOP();
	TIMING_PRINT("solvers_run()", "grid_clear");
	TIMING_START();

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic, 256)
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real pos_world[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		real pos_local[3] = {
			(pos_world[0] - bounds_min[0]) / bounds_wid[0],
			(pos_world[1] - bounds_min[1]) / bounds_wid[1],
			(pos_world[2] - bounds_min[2]) / bounds_wid[2]
		};

		for(unsigned int j = 0u; j < (unsigned int)GRID_SAMPLES; j++) {
			u32 s[4] = {
				(u32)0xB79ABC95u + (u32)j,
				(u32)0xE0DA3F84u + (u32)i,
				(u32)0xAB75F07Bu + (u32)step_number,
				(u32)0xCF52CA01u
			};

			pcg4d(s);
			pcg4d(s); // second round for better statistical quality

			real r[4] = {
				(real)INV_PCG32_MAX * (real)s[0],
				(real)INV_PCG32_MAX * (real)s[1],
				(real)INV_PCG32_MAX * (real)s[2],
				(real)INV_PCG32_MAX * (real)s[3]
			};

			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			real n[4] = {
				real_sqrt( (real)(-2.0) * real_log(r[0]) ) * real_cos( (real)TAU * r[1] ),
				real_sqrt( (real)(-2.0) * real_log(r[0]) ) * real_sin( (real)TAU * r[1] ),
				real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_cos( (real)TAU * r[3] ),
				real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_sin( (real)TAU * r[3] )
			};

			int coord[3] = {
				(int)( (real)0.500 * n[0] + (real)grid_dim[0] * pos_local[0] + (real)0.5),
				(int)( (real)0.500 * n[1] + (real)grid_dim[1] * pos_local[1] + (real)0.5),
				(int)( (real)0.500 * n[2] + (real)grid_dim[2] * pos_local[2] + (real)0.5)
			};

			if(0 <= coord[0] && coord[0] < (int)grid_dim[0]
			&& 0 <= coord[1] && coord[1] < (int)grid_dim[1]
			&& 0 <= coord[2] && coord[2] < (int)grid_dim[2]) {
				unsigned int index = grid_dim[1]*grid_dim[0]*(unsigned int)coord[2]+grid_dim[0]*(unsigned int)coord[1]+(unsigned int)coord[0];

				#ifdef _OPENMP
				#pragma omp atomic update
				#endif
				grid_accum[index] += (i32)1;
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("solvers_run()", "grid_dist");
	TIMING_START();

	double inv_grid_samples = 1.0 / (double)GRID_SAMPLES;

	for(unsigned int i = 0u; i < grid.dim[2u] * grid.dim[1u] * grid.dim[0u]; i++) {
		i32 val = (i32)0;

		#ifdef _OPENMP
		#pragma omp atomic read
		#endif
		val = grid_accum[i];

		grid_poten[i] = (f32)(inv_grid_samples * (double)val);
	}

	TIMING_STOP();
	TIMING_PRINT("solvers_run()", "grid_read");
	TIMING_START();

	grid_solve_fft(&grid);

	TIMING_STOP();
	TIMING_PRINT("solvers_run()", "solve_fft");
	#endif

	#ifdef PAIRWISE_SOLVER_DECIMATION
	unsigned int j_length = (          N / (unsigned int)PAIRWISE_SOLVER_DECIMATION);
	unsigned int j_offset = (step_number % (unsigned int)PAIRWISE_SOLVER_DECIMATION) * j_length;
	#endif

	TIMING_START();

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_SIM
		#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N], pot[:N]) map(to: pos[:3u*N], mas[:N])
		#else
		#pragma omp parallel for schedule(dynamic, 256)
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
		#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
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

			// body j's position
			real r_j[3] = {
				pos[3u*j+0u],
				pos[3u*j+1u],
				pos[3u*j+2u]
			};

			// vector to body i's position from body j's position
			real r_ij[3] = {
				r_i[0] - r_j[0],
				r_i[1] - r_j[1],
				r_i[2] - r_j[2]
			};

			// squared distance between bodies i and j
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
	}

	TIMING_STOP();
	TIMING_PRINT("solvers_run()", "brute_force");

	return STARFLOOD_SUCCESS;
}
