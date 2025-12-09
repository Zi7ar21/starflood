#include "simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "config.h"
#include "rng.h"

int simulation_init(simulation_t* simulation, unsigned int N) {
	simulation_t sim = *simulation;

	printf("Simulation Memory Addresses:\n");

	size_t pot_size = (size_t)1u * (size_t)N;
	size_t kin_size = (size_t)1u * (size_t)N;
	size_t mas_size = (size_t)1u * (size_t)N;
	size_t pos_size = (size_t)3u * (size_t)N;
	size_t vel_size = (size_t)3u * (size_t)N;
	size_t acc_size = (size_t)3u * (size_t)N;

	size_t pot_offset = (size_t)0u;
	size_t kin_offset = pot_offset + pot_size;
	size_t mas_offset = kin_offset + kin_size;
	size_t pos_offset = mas_offset + mas_size;
	size_t vel_offset = pos_offset + pos_size;
	size_t acc_offset = vel_offset + vel_size;

	size_t buf_size = acc_offset + acc_size;

	size_t mem_size = sizeof(real) * buf_size;

	void* mem = aligned_alloc(STARFLOOD_ALIGNMENT, mem_size);

	if(NULL == mem) {
		fprintf(stderr, "error in aligned_alloc(%zu, %zu) while allocating memory for the simulation", STARFLOOD_ALIGNMENT, mem_size);

		perror("");

		return EXIT_FAILURE;
	}

	printf("  mem: %p\n", mem);

	memset(mem, 0, mem_size);

	real* buf = (real*)mem;

	real* pot = &(buf[pot_offset]);
	real* kin = &(buf[kin_offset]);
	real* mas = &(buf[mas_offset]);
	real* pos = &(buf[pos_offset]);
	real* vel = &(buf[vel_offset]);
	real* acc = &(buf[acc_offset]);

	printf("  pot: %p (+%zu)\n", (void*)pot, pot_offset);
	printf("  kin: %p (+%zu)\n", (void*)kin, kin_offset);
	printf("  mas: %p (+%zu)\n", (void*)mas, mas_offset);
	printf("  pos: %p (+%zu)\n", (void*)pos, pos_offset);
	printf("  vel: %p (+%zu)\n", (void*)vel, vel_offset);
	printf("  acc: %p (+%zu)\n", (void*)acc, acc_offset);
	printf("\n");

	for(unsigned int i = 0u; i < N; i++) {
		pot[i] = (real)0.0;
	}

	for(unsigned int i = 0u; i < N; i++) {
		kin[i] = (real)0.0;
	}

	{
		real body_mass = (real)(10.0 / (double)N);

		for(unsigned int i = 0u; i < N; i++) {
			mas[i] = body_mass;
		}
	}

	for(unsigned int i = (size_t)0u; i < N; i++) {
		uint32_t s[4] = {(uint32_t)i, (uint32_t)420u, (uint32_t)69u, (uint32_t)1337u};

		pcg4d(s);
		pcg4d(s); // second round for better statistical quality

		{
			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			double p[3] = {
				sqrt( -2.0 * log(r[0]) ) * cos(TAU * r[2]),
				sqrt( -2.0 * log(r[0]) ) * sin(TAU * r[2]),
				sqrt( -2.0 * log(r[1]) ) * cos(TAU * r[3])
			};

			/*
			pos[3u*i+0u] = 2.0*((double)s[0]/(double)0xFFFFFFFFu)-1.0;
			pos[3u*i+1u] = 2.0*((double)s[1]/(double)0xFFFFFFFFu)-1.0;
			pos[3u*i+2u] = 2.0*((double)s[2]/(double)0xFFFFFFFFu)-1.0;
			*/

			pos[3u*i+0u] = (real)(1.000 * p[0]);
			pos[3u*i+1u] = (real)(0.100 * p[1]);
			pos[3u*i+2u] = (real)(1.000 * p[2]);
		}

		pcg4d(s);

		{
			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			double v[3] = {
				sqrt( -2.0 * log(r[0]) ) * cos(TAU * r[2]),
				sqrt( -2.0 * log(r[0]) ) * sin(TAU * r[2]),
				sqrt( -2.0 * log(r[1]) ) * cos(TAU * r[3])
			};

			vel[3u*i+0u] = (real)( INV_TAU * 0.200 * pos[3u*i+2u] + 0.0001 * v[0]);
			vel[3u*i+1u] = (real)( INV_TAU * 0.100 * pos[3u*i+1u] + 0.0001 * v[1]);
			vel[3u*i+2u] = (real)(-INV_TAU * 0.200 * pos[3u*i+0u] + 0.0001 * v[2]);
		}
	}

	for(size_t i = (size_t)0u; i < (size_t)3u*N; i++) {
		acc[i] = (real)0.0;
	}

	sim.step_number = 0u;

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

	const real dt = (real)TIMESTEP_SIZE;

	unsigned int N = sim.N;
	unsigned int step_number = sim.step_number;
	real* pot = sim.pot;
	real* kin = sim.kin;
	real* mas = sim.mas;
	real* pos = sim.pos;
	real* vel = sim.vel;
	real* acc = sim.acc;

	#ifdef _OPENMP
	//#pragma omp parallel for schedule(dynamic, 128)
	#pragma omp target teams distribute parallel for map(pot[0u:N], kin[0u:N], mas[0u:N], pos[0u:3u*N], vel[0u:3u*N], acc[0u:3u*N])
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real U_sum = (real)0.0;

		real U_c = (real)0.0;

		real F_sum[3] = {
			(real)0.0,
			(real)0.0,
			(real)0.0
		};

		real F_c[3] = {
			(real)0.0,
			(real)0.0,
			(real)0.0
		};

		real m_i = mas[i]; // Body mass

		real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		for(unsigned int j = 0u; j < N; j++) {
		//for(unsigned int j = (N/4u)*(step_number % 4u); j < (N/4u)*(step_number % 4u) + (N/4u); j++) {
			if(i == j) {
				continue;
			}

			real m_j = mas[j];

			real r_j[3] = {
				pos[3u*j+0u],
				pos[3u*j+1u],
				pos[3u*j+2u]
			};

			real r_ij[3] = {
				r_j[0u] - r_i[0u],
				r_j[1u] - r_i[1u],
				r_j[2u] - r_i[2u]
			};

			real r2 = (r_ij[0u]*r_ij[0u])+(r_ij[1u]*r_ij[1u])+(r_ij[2u]*r_ij[2u]);

			real inv_r2 = (real)1.0 / (      r2 +(real)0.000001);
			real inv_r1 = (real)1.0 / (sqrtf(r2)+(real)0.000001);

			real U_ij = -(real)G * m_i * m_j * inv_r2;

			{
				// Naïve summation
				//U_sum += U_ij;

				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = U_ij - U_c;
				volatile real t = U_sum + y;
				volatile real z = t - U_sum;
				U_c = z - y;
				U_sum = t;
			}

			real F_ij[3] = {
				-U_ij * r_ij[0u] * inv_r1,
				-U_ij * r_ij[1u] * inv_r1,
				-U_ij * r_ij[2u] * inv_r1
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				// Naïve summation
				//F_sum[k] += F_ij[k];

				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				real y = F_ij[k] - F_c[k];
				volatile real t = F_sum[k] + y;
				volatile real z = t - F_sum[k];
				F_c[k] = z - y;
				F_sum[k] = t;
			}
		}

		pot[i] = U_sum;

		acc[3u*i+0u] = F_sum[0u];
		acc[3u*i+1u] = F_sum[1u];
		acc[3u*i+2u] = F_sum[2u];
	}

	for(unsigned int i = 0u; i < N; i++) {
		real v_i[3] = {
			(real)vel[3u*i+0u],
			(real)vel[3u*i+1u],
			(real)vel[3u*i+2u]
		};

		kin[i] = (real)0.5 * mas[i] * ((v_i[0]*v_i[0])+(v_i[1]*v_i[1])+(v_i[2]*v_i[2])); // K = (1/2) * m * v^2
	}

	// Compute total energy
	{
		double U_sum = 0.0, U_c = 0.0;
		double K_sum = 0.0, K_c = 0.0;

		for(unsigned int i = 0u; i < N; i++) {
			// Naïve summation
			//U_sum += (double)pot[i];

			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			double y = (double)pot[i] - U_c;
			volatile double t = U_sum + y;
			volatile double z = t - U_sum;
			U_c = z - y;
			U_sum = t;
		}

		for(unsigned int i = 0u; i < N; i++) {
			// Naïve summation
			//K_sum += (double)kin[i];

			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			double y = (double)kin[i] - K_c;
			volatile double t = K_sum + y;
			volatile double z = t - K_sum;
			K_c = z - y;
			K_sum = t;
		}

		printf("E (E_tot) = % .015f\nU (E_pot) = % .015f\nK (E_kin) = % .015f\n", U_sum + K_sum, U_sum, K_sum);
	}

	// 1/2 kick
	for(unsigned int i = 0u; i < 3u * N; i++) {
		pos[i] += (real)0.5 * dt * vel[i];
	}

	// Drift
	for(unsigned int i = 0u; i < 3u * N; i++) {
		vel[i] += dt * acc[i];
	}

	// 1/2 kick
	for(unsigned int i = 0u; i < 3u * N; i++) {
		pos[i] += (real)0.5 * dt * vel[i];
	}

	sim.step_number++;

	*simulation = sim;

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

int simulation_read(simulation_t* simulation, unsigned int N, const char* restrict filename) {
	simulation_t sim = *simulation;

	if( sizeof(real) != sizeof(float) ) {
		return EXIT_FAILURE;
	}

	FILE* file = fopen(filename, "rb");

	if(NULL == file) {
		fprintf(stderr, "Error: fopen(\"%s\", \"%s\") ", filename, "rb");

		perror("failed");

		return EXIT_FAILURE;
	}

	/*
	fseek(file, (long)2, SEEK_SET);

	const uint32_t starflood_magic_number = (uint32_t)0x40C90FDBu;

	uint32_t magic_number = (uint32_t)0xDEADBEEFu;

	fread(&magic_number, sizeof(uint32_t), (size_t)1u, file);

	if(starflood_magic_number!= magic_number) {
		fprintf(stderr, "Error: Expected magic number 0x%08X, got 0x%08X instead!\n", starflood_magic_number, magic_number);

		fclose(file);

		return EXIT_FAILURE;
	}

	uint32_t N = (uint32_t)0u;

	fread(&N, sizeof(uint32_t), (size_t)1u, file);

	if((uint32_t)(1024u*1024u) < N) {
		fprintf(stderr, "Error: Assuming N = %u is too large!\n", (unsigned int)N);

		fclose(file);

		return EXIT_FAILURE;
	}

	fseek(file, (long)10, SEEK_SET);
	*/

	if( EXIT_SUCCESS != simulation_init(&sim, N) ) {
		fclose(file);

		return EXIT_FAILURE;
	}

	//fread(sim.mem, sizeof(real), (size_t)N * (size_t)(3u*1u+3u*3u), file);

	fread(sim.mas, sizeof(float), (size_t)N * (size_t)(1u), file);
	fread(sim.pos, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.vel, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.acc, sizeof(float), (size_t)N * (size_t)(3u), file);

	fclose(file);

	*simulation = sim;

	return EXIT_SUCCESS;
}

int simulation_dump(simulation_t* simulation, const char* restrict filename) {
	simulation_t sim = *simulation;

	FILE* file = fopen(filename, "wb");

	if(NULL == file) {
		fprintf(stderr, "Error: fopen(\"%s\", \"%s\") ", filename, "wb");

		perror("failed");

		return EXIT_FAILURE;
	};

	uint32_t N = (uint32_t)sim.N;

	/*
	const uint8_t magic_bytes[2] = {(uint8_t)0x53u, (uint8_t)0x46u};

	fwrite(magic_bytes, sizeof(uint8_t), (size_t)2u, file);

	const uint32_t starflood_magic_number = (uint32_t)0x40C90FDBu;

	fwrite(&starflood_magic_number, sizeof(uint32_t), (size_t)1u, file);

	fwrite(&N, sizeof(uint32_t), (size_t)1u, file);
	*/

	fwrite(sim.mem, sizeof(real), (size_t)N * (size_t)(3u*1u+3u*3u), file);

	fclose(file);

	return EXIT_SUCCESS;
}
