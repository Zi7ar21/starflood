// needed for clock_gettime()
#define _POSIX_C_SOURCE 200112L

#include "simulation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "config.h"
#include "initcond.h"
#include "log.h"
#include "rng.h"
#include "solver.h"
#include "timing.h"

#ifdef LOG_STATISTICS
log_t log_statistics;
#endif

#ifdef LOG_TIMINGS_SIM_STEP
log_t log_timings_sim_step;
#endif

int simulation_init(simulation_t* restrict simulation, unsigned int N) {
	TIMING_INIT();

	#ifdef LOG_STATISTICS
	if( STARFLOOD_SUCCESS != log_init(&log_statistics, LOG_STATISTICS) ) {
		fprintf(stderr, "%s error: %s failed!\n", "simulation_init()", "log_init(&log_statistics)");
		return STARFLOOD_FAILURE;
	}

	fprintf(log_statistics.file, "%s,%s,%s,%s\n", "step_number", "total_energy", "potential_energy", "kinetic_energy");
	fflush(log_statistics.file);
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	if( STARFLOOD_SUCCESS != log_init(&log_timings_sim_step, LOG_TIMINGS_SIM_STEP) ) {
		fprintf(stderr, "%s error: %s failed!\n", "simulation_init()", "log_init(&log_timings_sim_step)");
		return STARFLOOD_FAILURE;
	}

	fprintf(log_timings_sim_step.file, "%s,%s,%s,%s,%s,%s,%s,%s\n", "step_number", "kick_0", "drift", "solver_run", "calc_kin", "pot_sum", "kin_sum", "kick_1");
	fflush(log_timings_sim_step.file);
	#endif

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

	TIMING_START();

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
	#endif

	TIMING_STOP();

	#ifdef STARFLOOD_ALIGNMENT
	TIMING_PRINT("simulation_init()", "posix_memalign()");
	#else
	TIMING_PRINT("simulation_init()", "malloc()");
	#endif

	if(NULL == mem) {
		#ifdef STARFLOOD_ALIGNMENT
		fprintf(stderr, "%s error: mem is NULL after posix_memalign(&mem, %zu, %zu", "simulation_init()", (size_t)STARFLOOD_ALIGNMENT, mem_size);
		#else
		fprintf(stderr, "%s error: mem is NULL after malloc(%zu", "simulation_init()", mem_size);
		#endif

		perror(")");

		return STARFLOOD_FAILURE;
	}

	printf("  mem: %p\n", mem);

	TIMING_START();

	memset(mem, 0, mem_size);

	TIMING_STOP();
	TIMING_PRINT("simulation_init()", "memset()");

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

	TIMING_START();

	#ifdef ENABLE_SIMULATION
	initcond_generate(mas, pos, vel, N);
	#endif

	TIMING_STOP();
	TIMING_PRINT("simulation_init()", "initcond_generate");
	TIMING_START();

	#ifdef ENABLE_SIMULATION
	solver_run(pot, acc, mas, pos, N, 0u);
	#endif

	TIMING_STOP();
	TIMING_PRINT("simulation_init()", "solver_run");
	TIMING_START();

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

	return STARFLOOD_SUCCESS;
}

int simulation_free(simulation_t* restrict simulation) {
	TIMING_INIT();

	#ifdef LOG_STATISTICS
	if( STARFLOOD_SUCCESS != log_free(&log_statistics) ) {
		fprintf(stderr, "%s error: %s failed!\n", "simulation_free()", "log_free(&log_statistics)");
	}
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	if( STARFLOOD_SUCCESS != log_free(&log_timings_sim_step) ) {
		fprintf(stderr, "%s error: %s failed!\n", "simulation_free()", "log_free(&log_timings_sim_step)");
	}
	#endif

	simulation_t sim = *simulation;

	TIMING_START();

	free(sim.mem);

	TIMING_STOP();
	TIMING_PRINT("simulation_free()", "free()");

	sim.mem = NULL;
	sim.pot = (real*)NULL;
	sim.kin = (real*)NULL;
	sim.mas = (real*)NULL;
	sim.pos = (real*)NULL;
	sim.vel = (real*)NULL;
	sim.acc = (real*)NULL;

	*simulation = sim;

	return STARFLOOD_SUCCESS;
}

int simulation_read(simulation_t* restrict simulation, const char* restrict filename) {
	TIMING_INIT();

	simulation_t sim = *simulation;

	if( sizeof(real) != sizeof(float) ) {
		return STARFLOOD_FAILURE;
	}

	unsigned int N = sim.N;

	if( 0u >= N ) {
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "rb");

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fopen()");

	if(NULL == file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "simulation_read()", filename, "rb");

		perror("failed");

		return STARFLOOD_FAILURE;
	}

	/*
	if( STARFLOOD_SUCCESS != simulation_init(&sim, N) ) {
		fclose(file);

		return STARFLOOD_FAILURE;
	}
	*/

	TIMING_START();

	fread(sim.mem, sizeof(real), (size_t)N * (size_t)(3u*1u+3u*3u), file);

	/*
	fread(sim.mas, sizeof(float), (size_t)N * (size_t)(1u), file);
	fread(sim.pos, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.vel, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.acc, sizeof(float), (size_t)N * (size_t)(3u), file);
	*/

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fread()");
	TIMING_START();

	fclose(file);

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fclose()");

	*simulation = sim;

	return STARFLOOD_SUCCESS;
}

int simulation_save(simulation_t* restrict simulation, const char* restrict filename) {
	TIMING_INIT();

	simulation_t sim = *simulation;

	uint32_t N = (uint32_t)sim.N;

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fopen()");
	TIMING_START();

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "simulation_save()", filename, "wb");

		perror("failed");

		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	/*
	const uint8_t magic_bytes[2] = {(uint8_t)0x53u, (uint8_t)0x46u};

	fwrite(magic_bytes, sizeof(uint8_t), (size_t)2u, file);

	const uint32_t starflood_magic_number = (uint32_t)0x40C90FDBu;

	fwrite(&starflood_magic_number, sizeof(uint32_t), (size_t)1u, file);

	fwrite(&N, sizeof(uint32_t), (size_t)1u, file);
	*/

	fwrite(sim.mem, sizeof(real), (size_t)N * (size_t)(3u*1u+3u*3u), file);

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fwrite()");
	TIMING_START();

	if( 0 != fclose(file) ) {
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fclose()");
	TIMING_START();

	return STARFLOOD_SUCCESS;
}

int simulation_step(simulation_t* restrict simulation) {
	TIMING_INIT();

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

	// for the very first leapfrog kick
	//if(0u >= step_number) {
	//	solver_run(pot, acc, mas, pos, N, step_number);
	//}

	#ifdef LOG_STATISTICS
	fprintf(log_statistics.file, "%u", step_number);
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	fprintf(log_timings_sim_step.file, "%u", step_number);
	#endif

	TIMING_START();

	// kick_0
	for(unsigned int i = 0u; i < 3u * N; i++) {
		vel[i] += (real)0.5 * dt * acc[i];
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_step()", "kick_0");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif
	TIMING_START();

	// drift
	for(unsigned int i = 0u; i < 3u * N; i++) {
		pos[i] += dt * vel[i];
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_step()", "drift");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif
	TIMING_START();

	// run solver
	solver_run(pot, acc, mas, pos, N, step_number);

	TIMING_STOP();
	TIMING_PRINT("simulation_step()", "solver_run");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif
	TIMING_START();

	// recalculate kinetic energy
	for(unsigned int i = 0u; i < N; i++) {
		real m_i = mas[i];

		real v_i[3] = {
			vel[3u*i+0u],
			vel[3u*i+1u],
			vel[3u*i+2u]
		};

		kin[i] = (real)0.5 * m_i * ((v_i[0]*v_i[0])+(v_i[1]*v_i[1])+(v_i[2]*v_i[2])); // K = (1/2) * m * v^2
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_step()", "calc_kin");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif

	// compute total energy
	{
		double U_sum = 0.0, U_c = 0.0;
		double K_sum = 0.0, K_c = 0.0;

		TIMING_START();

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

		TIMING_STOP();
		TIMING_PRINT("simulation_step()", "pot_sum");
		#ifdef LOG_TIMINGS_SIM_STEP
		LOG_TIMING(log_timings_sim_step);
		#endif
		TIMING_START();

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

		TIMING_STOP();
		TIMING_PRINT("simulation_step()", "kin_sum");
		#ifdef LOG_TIMINGS_SIM_STEP
		LOG_TIMING(log_timings_sim_step);
		#endif

		#ifdef LOG_STATISTICS
		fprintf(log_statistics.file, ",%.015f,%.015f,%.015f", U_sum + K_sum, U_sum, K_sum);
		#endif

		printf("E (E_tot) = % .015f\nU (E_pot) = % .015f\nK (E_kin) = % .015f\n", U_sum + K_sum, U_sum, K_sum);
	}

	TIMING_START();

	// kick_1
	for(unsigned int i = 0u; i < 3u * N; i++) {
		vel[i] += (real)0.5 * dt * acc[i];
	}

	TIMING_STOP();	
	TIMING_PRINT("simulation_step()", "kick_1");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif

	#ifdef LOG_STATISTICS
	fprintf(log_statistics.file, "\n");
	fflush(log_statistics.file);
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	fprintf(log_timings_sim_step.file, "\n");
	fflush(log_timings_sim_step.file);
	#endif

	*simulation = sim;

	return STARFLOOD_SUCCESS;
}
