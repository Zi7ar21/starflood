// Needed for posix_memalign()
#define _POSIX_C_SOURCE 200112L

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
#include "solver.h"

int simulation_init(simulation_t* simulation, unsigned int N) {
	simulation_t sim = *simulation;

	void* mem = NULL;

	real* pot = (real*)NULL;
	real* kin = (real*)NULL;
	real* mas = (real*)NULL;
	real* pos = (real*)NULL;
	real* vel = (real*)NULL;
	real* acc = (real*)NULL;

	printf("Simulation Memory Addresses:\n");

	size_t pot_length = (size_t)1u * (size_t)N;
	size_t kin_length = (size_t)1u * (size_t)N;
	size_t mas_length = (size_t)1u * (size_t)N;
	size_t pos_length = (size_t)3u * (size_t)N;
	size_t vel_length = (size_t)3u * (size_t)N;
	size_t acc_length = (size_t)3u * (size_t)N;

	size_t pot_offset = (size_t)0u;
	size_t kin_offset = pot_offset + pot_length;
	size_t mas_offset = kin_offset + kin_length;
	size_t pos_offset = mas_offset + mas_length;
	size_t vel_offset = pos_offset + pos_length;
	size_t acc_offset = vel_offset + vel_length;

	size_t mem_size = sizeof(real) * (acc_offset + acc_length);

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
	#endif

	if(NULL == mem) {
		#ifdef STARFLOOD_ALIGNMENT
		fprintf(stderr, "error in posix_memalign(&mem, %zu, %zu) while allocating memory for the ", (size_t)STARFLOOD_ALIGNMENT, mem_size);
		#else
		fprintf(stderr, "error in malloc(%zu) while allocating memory for the ", mem_size);
		#endif

		perror("simulation");

		return EXIT_FAILURE;
	}

	printf("  mem: %p\n", mem);

	memset(mem, 0, mem_size);

	pot = (real*)mem + (size_t)pot_offset;
	kin = (real*)mem + (size_t)kin_offset;
	mas = (real*)mem + (size_t)mas_offset;
	pos = (real*)mem + (size_t)pos_offset;
	vel = (real*)mem + (size_t)vel_offset;
	acc = (real*)mem + (size_t)acc_offset;

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

	double body_mass = 1.0 / (double)N;

	for(unsigned int i = 0u; i < N; i++) {
		mas[i] = (real)body_mass;
	}

	for(unsigned int i = 0u; i < N; i++) {
		double p[3] = {0.0, 0.0, 0.0}; // position
		double v[3] = {0.0, 0.0, 0.0}; // velocity

		// Initialize position
		{
			uint32_t s[4] = {
				(uint32_t)i,
				(uint32_t)0x203F83D3u,
				(uint32_t)0x710BE493u,
				(uint32_t)0xA294E7F1u
			};

			pcg4d(s);
			pcg4d(s); // second round for better statistical quality

			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			double n[3] = {
				sqrt( -2.0 * log(r[0u]) ) * cos(TAU * r[2u]),
				sqrt( -2.0 * log(r[0u]) ) * sin(TAU * r[2u]),
				sqrt( -2.0 * log(r[1u]) ) * cos(TAU * r[3u])
			};

			/*
			p[0u] = 1.000 * n[0u];
			p[1u] = 0.100 * n[1u];
			p[2u] = 1.000 * n[2u];
			*/

			p[0u] = 1.000 * n[0u];
			p[1u] = 1.000 * n[1u];
			p[2u] = 1.000 * n[2u];
		}

		// Initialize velocity
		{
			uint32_t s[4] = {
				(uint32_t)i,
				(uint32_t)0xAE8D7CF2u,
				(uint32_t)0xE4827F00u,
				(uint32_t)0xE44CA389u
			};

			pcg4d(s);
			pcg4d(s); // second round for better statistical quality

			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			double n[3] = {
				sqrt( -2.0 * log(r[0u]) ) * cos(TAU * r[2u]),
				sqrt( -2.0 * log(r[0u]) ) * sin(TAU * r[2u]),
				sqrt( -2.0 * log(r[1u]) ) * cos(TAU * r[3u])
			};

			/*
			v[0u] =  INV_TAU * 0.125 * p[2u] + 0.001 * n[0u];
			v[1u] =  INV_TAU * 0.000 * p[1u] + 0.001 * n[1u];
			v[2u] = -INV_TAU * 0.125 * p[0u] + 0.001 * n[2u];
			*/

			/*
			double r2 = (p[0u]*p[0u])+(p[2u]*p[2u]);

			double inv_r2 = 1.0 / (      r2+0.1  );
			double inv_r1 = 1.0 / ( sqrt(r2+0.1) );
			*/

			v[0u] = 0.250 * sqrt(G * body_mass) *  p[2u] + 0.000001 * n[0u];
			v[1u] = 0.000 * sqrt(G * body_mass) *  p[1u] + 0.000001 * n[1u];
			v[2u] = 0.250 * sqrt(G * body_mass) * -p[0u] + 0.000001 * n[2u];
		}

		{
			uint32_t s[4] = {
				(uint32_t)i,
				(uint32_t)0xE5C4E174u,
				(uint32_t)0xCB77B51Eu,
				(uint32_t)0xA82C87A3u
			};

			pcg4d(s);
			pcg4d(s); // second round for better statistical quality

			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Thin disk
			if(0.000 <= r[0u] && r[0u] < 0.850) {
				p[1u] *= 0.025;
			}

			/*
			p[0u] += r[1u] < 0.500 ?  4.000 : -4.000;
			p[1u] += r[1u] < 0.500 ?  1.000 : -1.000;
			v[0u] *= r[1u] < 0.500 ?  1.000 : -1.000;
			v[2u] *= r[1u] < 0.500 ?  1.000 : -1.000;
			v[0u] += r[1u] < 0.500 ? -0.005 :  0.005;
			v[1u] += r[1u] < 0.500 ? -0.005 :  0.005;
			*/
		}

		pos[3u*i+0u] = (real)p[0u];
		pos[3u*i+1u] = (real)p[1u];
		pos[3u*i+2u] = (real)p[2u];

		vel[3u*i+0u] = (real)v[0u];
		vel[3u*i+1u] = (real)v[1u];
		vel[3u*i+2u] = (real)v[2u];
	}

	for(unsigned int i = 0u; i < 3u * N; i++) {
		acc[i] = (real)0.0;
	}

	#ifdef SIMULATION_ENABLED
	solver_run(pot, acc, mas, pos, N, 0u);
	#endif

	sim.N   = N;
	sim.step_number = 0u;
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
	#ifdef _OPENMP
	double t0 = omp_get_wtime();
	double t1 = omp_get_wtime();
	#endif

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
	t0 = omp_get_wtime();
	#endif

	// Kick
	for(unsigned int i = 0u; i < 3u * N; i++) {
		vel[i] += (real)0.5 * dt * acc[i];
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("simulation_step: %.09f ms kick\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	// Drift
	for(unsigned int i = 0u; i < 3u * N; i++) {
		pos[i] += dt * vel[i];
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("simulation_step: %.09f ms drift\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	// Update potential energy/acceleration
	solver_run(pot, acc, mas, pos, N, step_number);

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("simulation_step: %.09f ms solver_run\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	// Update kinetic energy
	for(unsigned int i = 0u; i < N; i++) {
		real v_i[3] = {
			(real)vel[3u*i+0u],
			(real)vel[3u*i+1u],
			(real)vel[3u*i+2u]
		};

		kin[i] = (real)0.5 * mas[i] * ((v_i[0]*v_i[0])+(v_i[1]*v_i[1])+(v_i[2]*v_i[2])); // K = (1/2) * m * v^2
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("simulation_step: %.09f ms update kinetic energy\n", 1000.0*(t1-t0));
	#endif

	// Compute total energy
	{
		#ifdef _OPENMP
		t0 = omp_get_wtime();
		#endif

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

		#ifdef _OPENMP
		t1 = omp_get_wtime();
		#endif

		printf("E (E_tot) = % .015f\nU (E_pot) = % .015f\nK (E_kin) = % .015f\n", U_sum + K_sum, U_sum, K_sum);
	}

	#ifdef _OPENMP
	printf("simulation_step: %.09f ms calculate energy sums\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	// Kick
	for(unsigned int i = 0u; i < 3u * N; i++) {
		vel[i] += (real)0.5 * dt * acc[i];
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("simulation_step: %.09f ms kick\n", 1000.0*(t1-t0));
	#endif

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
	sim.pot = (real*)NULL;
	sim.kin = (real*)NULL;
	sim.mas = (real*)NULL;
	sim.pos = (real*)NULL;
	sim.vel = (real*)NULL;
	sim.acc = (real*)NULL;

	*simulation = sim;

	return EXIT_SUCCESS;
}

int simulation_read(simulation_t* simulation, const char* restrict filename) {
	simulation_t sim = *simulation;

	if( sizeof(real) != sizeof(float) ) {
		return EXIT_FAILURE;
	}

	unsigned int N = sim.N;

	if( 0u >= N ) {
		return EXIT_FAILURE;
	}

	FILE* file = fopen(filename, "rb");

	if(NULL == file) {
		fprintf(stderr, "Error: fopen(\"%s\", \"%s\") ", filename, "rb");

		perror("failed");

		return EXIT_FAILURE;
	}

	/*
	if( EXIT_SUCCESS != simulation_init(&sim, N) ) {
		fclose(file);

		return EXIT_FAILURE;
	}
	*/

	fread(sim.mem, sizeof(real), (size_t)N * (size_t)(3u*1u+3u*3u), file);

	/*
	fread(sim.mas, sizeof(float), (size_t)N * (size_t)(1u), file);
	fread(sim.pos, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.vel, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.acc, sizeof(float), (size_t)N * (size_t)(3u), file);
	*/

	fclose(file);

	*simulation = sim;

	return EXIT_SUCCESS;
}

int simulation_save(simulation_t* simulation, const char* restrict filename) {
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
