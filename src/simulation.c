#include "simulation.h"

#include "pcg4d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int simulation_init(simulation_t* simulation, unsigned int N) {
	simulation_t sim = *simulation;

	size_t pot_size = (size_t)1u * N;
	size_t kin_size = (size_t)1u * N;
	size_t mas_size = (size_t)1u * N;
	size_t pos_size = (size_t)3u * N;
	size_t vel_size = (size_t)3u * N;
	size_t acc_size = (size_t)3u * N;

	size_t pot_offset = (size_t)0u;
	size_t kin_offset = pot_offset + pot_size;
	size_t mas_offset = kin_offset + kin_size;
	size_t pos_offset = mas_offset + mas_size;
	size_t vel_offset = pos_offset + pos_size;
	size_t acc_offset = vel_offset + vel_size;

	size_t buf_size = acc_offset + acc_size;

	size_t mem_size = sizeof(real)*buf_size;

	const size_t align_size = (size_t)4096u; // Set to page size or something similar

	void* mem = aligned_alloc(align_size, mem_size);

	if(NULL == mem) {
		fprintf(stderr, "error in aligned_alloc(%zu, %zu) while allocating memory for the simulation", align_size, mem_size);

		perror("");

		return EXIT_FAILURE;
	}

	printf("Shared Memory Address: %p\n", mem);

	memset(mem, 0, mem_size);

	real* buf = (real*)mem;

	real* pot = &(buf[pot_offset]);
	real* kin = &(buf[kin_offset]);
	real* mas = &(buf[mas_offset]);
	real* pos = &(buf[pos_offset]);
	real* vel = &(buf[vel_offset]);
	real* acc = &(buf[acc_offset]);

	printf("%p\n", (void*)pot);
	printf("%p\n", (void*)kin);
	printf("%p\n", (void*)mas);
	printf("%p\n", (void*)pos);
	printf("%p\n", (void*)vel);
	printf("%p\n", (void*)acc);

	for(size_t i = (size_t)0u; i < N; i++) {
		pot[i] = (real)0.0;
	}

	for(size_t i = (size_t)0u; i < N; i++) {
		kin[i] = (real)0.0;
	}

	for(size_t i = (size_t)0u; i < N; i++) {
		mas[i] = (double)1.0/(double)N;
	}

	for(size_t i = (size_t)0u; i < N; i++) {
		uint32_t s[4] = {(uint32_t)i, (uint32_t)420u, (uint32_t)69u, (uint32_t)1337u};

		pcg4d(s);

		pos[3u*i+0u] = 2.0*((double)s[0]/(double)0xFFFFFFFFu)-1.0;
		pos[3u*i+1u] = 2.0*((double)s[1]/(double)0xFFFFFFFFu)-1.0;
		pos[3u*i+2u] = 2.0*((double)s[2]/(double)0xFFFFFFFFu)-1.0;
	}

	for(size_t i = (size_t)0u; i < (size_t)3u*N; i++) {
		//vel[i] = (real)0.0;
		vel[i] = (real)0.001*pos[i];
	}

	for(size_t i = (size_t)0u; i < (size_t)3u*N; i++) {
		acc[i] = (real)0.0;
	}

	sim.N   = N;
	sim.mem = mem;
	sim.pot = pot;
	sim.kin = kin;
	sim.mas = mas;
	sim.pos = pos;
	sim.vel = vel;
	sim.acc = acc;

	*simulation = sim;

	return EXIT_SUCCESS;
}

int simulation_step(simulation_t* simulation) {
	simulation_t sim = *simulation;

	const real timestep = (real)0.3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068;
	const real G = (real)5.0;

	unsigned int N = sim.N;
	real* pot = sim.pot;
	real* kin = sim.kin;
	real* mas = sim.mas;
	real* pos = sim.pos;
	real* vel = sim.vel;
	real* acc = sim.acc;

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,100)
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real U_sum = (real)0;

		real F[3] = {(real)0.0, (real)0.0, (real)0.0};

		for(unsigned int j = 0u; j < N; j++) {
			if(i == j) {
				continue;
			}

			real dr[3] = {
				pos[3u*i+0u]-pos[3u*j+0u],
				pos[3u*i+1u]-pos[3u*j+1u],
				pos[3u*i+2u]-pos[3u*j+2u]
			};

			real r2 = (dr[0]*dr[0])+(dr[1]*dr[1])+(dr[2]*dr[2]);

			real r1 = (real)sqrtf(r2);

			real rinv = (real)1.0 / r1;

			real U = -G*mas[i]*mas[j]*rinv;

			U_sum += U;

			F[0] += U*dr[0]/(r2+(real)0.001);
			F[1] += U*dr[1]/(r2+(real)0.001);
			F[2] += U*dr[2]/(r2+(real)0.001);
		}

		pot[i] = U_sum;
		acc[3u*i+0u] = F[0];
		acc[3u*i+1u] = F[1];
		acc[3u*i+2u] = F[2];
	}

	for(unsigned int i = 0u; i < (size_t)N; i++) {
		kin[i] = (real)0.5*mas[i]*vel[i]*vel[i];
	}

	{
		double usum = 0.0;
		double ksum = 0.0;

		for(unsigned int i = 0u; i < N; i++) {
			usum += pot[i];
		}

		for(unsigned int i = 0u; i < N; i++) {
			ksum += kin[i];
		}

		printf("T = %.015f\nU = %.015f\nV = %.015f\n", usum+ksum, usum, ksum);
	}

	for(unsigned int i = 0u; i < (size_t)3u*N; i++) {
		pos[i] += (real)0.5*timestep*vel[i];
	}

	for(unsigned int i = 0u; i < (size_t)3u*N; i++) {
		vel[i] += timestep*acc[i];
	}

	for(unsigned int i = 0u; i < (size_t)3u*N; i++) {
		pos[i] += (real)0.5*timestep*vel[i];
	}

	//*simulation = sim;

	return EXIT_SUCCESS;
}

int simulation_free(simulation_t* simulation) {
	simulation_t sim = *simulation;

	{
		void* mem = sim.mem;

		free(mem);
	}

	sim.mem = NULL;
	sim.pot = NULL;
	sim.kin = NULL;
	sim.mas = NULL;
	sim.pos = NULL;
	sim.vel = NULL;
	sim.acc = NULL;

	*simulation = sim;

	return EXIT_SUCCESS;
}

