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

int simulation_init(sim_t* restrict simulation, unsigned int N) {
	TIMING_INIT();

	#ifdef LOG_STATISTICS
	if( STARFLOOD_SUCCESS != log_init(&log_statistics, OUTPUT_DIR "/" LOG_STATISTICS) ) {
		fprintf(stderr, "%s error: %s failed!\n", "simulation_init()", "log_init(&log_statistics)");
		return STARFLOOD_FAILURE;
	}

	fprintf(log_statistics.file, "%s,%s,%s,%s\n", "step_number", "total_energy", "potential_energy", "kinetic_energy");
	fflush(log_statistics.file);
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	if( STARFLOOD_SUCCESS != log_init(&log_timings_sim_step, OUTPUT_DIR "/" LOG_TIMINGS_SIM_STEP) ) {
		fprintf(stderr, "%s error: %s failed!\n", "simulation_init()", "log_init(&log_timings_sim_step)");
		return STARFLOOD_FAILURE;
	}

	fprintf(log_timings_sim_step.file, "%s,%s,%s,%s,%s,%s,%s\n", "step_number", "kick_0", "drift", "solver_run", "ken_sum", "pen_sum", "kick_1");
	fflush(log_timings_sim_step.file);
	#endif

	sim_t sim = *simulation;

	void* mem = NULL;

	real* mas = (real*)NULL;
	real* rad = (real*)NULL;
	real* pos = (real*)NULL;
	real* vel = (real*)NULL;

	real* acc = (real*)NULL;

	real* pot = (real*)NULL;

	real* rho = (real*)NULL;
	real* prs = (real*)NULL;

	real* ken = (real*)NULL;
	real* pen = (real*)NULL;

	printf("Simulation Memory Addresses:\n");

	// calculate buffer lengths
	size_t mas_length = (size_t)1u * (size_t)N;
	#ifdef ENABLE_SPH
	size_t rad_length = (size_t)1u * (size_t)N;
	#endif
	size_t pos_length = (size_t)3u * (size_t)N;
	size_t vel_length = (size_t)3u * (size_t)N;
	size_t acc_length = (size_t)3u * (size_t)N;
	size_t pot_length = (size_t)1u * (size_t)N;
	#ifdef ENABLE_SPH
	size_t rho_length = (size_t)1u * (size_t)N;
	size_t prs_length = (size_t)1u * (size_t)N;
	#endif
	size_t ken_length = (size_t)1u * (size_t)N;
	size_t pen_length = (size_t)1u * (size_t)N;

	// calculate buffer offsets
	size_t mem_length = (size_t)0u;

	size_t mas_offset = mem_length;
	mem_length += mas_length;

	#ifdef ENABLE_SPH
	size_t rad_offset = mem_length;
	mem_length += rad_length;
	#endif

	size_t pos_offset = mem_length;
	mem_length += pos_length;

	size_t vel_offset = mem_length;
	mem_length += vel_length;

	size_t acc_offset = mem_length;
	mem_length += acc_length;

	size_t pot_offset = mem_length;
	mem_length += pot_length;

	#ifdef ENABLE_SPH
	size_t rho_offset = mem_length;
	mem_length += rho_length;

	size_t prs_offset = mem_length;
	mem_length += prs_length;
	#endif

	size_t ken_offset = mem_length;
	mem_length += ken_length;

	size_t pen_offset = mem_length;
	mem_length += pen_length;

	size_t mem_size = sizeof(real) * mem_length;

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

	printf("  mem: %p (%zu bytes)\n", mem, mem_size);

	TIMING_START();

	memset(mem, 0, mem_size);

	TIMING_STOP();
	TIMING_PRINT("simulation_init()", "memset()");

	{
		mas = &( (real*)mem )[mas_offset];
		#ifdef ENABLE_SPH
		rad = &( (real*)mem )[rad_offset];
		#endif
		pos = &( (real*)mem )[pos_offset];
		vel = &( (real*)mem )[vel_offset];

		acc = &( (real*)mem )[acc_offset];

		pot = &( (real*)mem )[pot_offset];

		#ifdef ENABLE_SPH
		rho = &( (real*)mem )[rho_offset];
		prs = &( (real*)mem )[prs_offset];
		#endif

		ken = &( (real*)mem )[ken_offset];
		pen = &( (real*)mem )[pen_offset];
	}

	printf("  mas: %p (%zu bytes)\n", (void*)mas, sizeof(real)*(size_t)mas_length);
	#ifdef ENABLE_SPH
	printf("  rad: %p (%zu bytes)\n", (void*)rad, sizeof(real)*(size_t)rad_length);
	#endif
	printf("  pos: %p (%zu bytes)\n", (void*)pos, sizeof(real)*(size_t)pos_length);
	printf("  vel: %p (%zu bytes)\n", (void*)vel, sizeof(real)*(size_t)vel_length);
	printf("  acc: %p (%zu bytes)\n", (void*)acc, sizeof(real)*(size_t)acc_length);
	printf("  pot: %p (%zu bytes)\n", (void*)pot, sizeof(real)*(size_t)pot_length);
	#ifdef ENABLE_SPH
	printf("  rho: %p (%zu bytes)\n", (void*)rho, sizeof(real)*(size_t)rho_length);
	printf("  prs: %p (%zu bytes)\n", (void*)prs, sizeof(real)*(size_t)prs_length);
	#endif
	printf("  ken: %p (%zu bytes)\n", (void*)ken, sizeof(real)*(size_t)ken_length);
	printf("  pen: %p (%zu bytes)\n", (void*)pen, sizeof(real)*(size_t)pen_length);
	printf("\n");

	TIMING_START();

	#ifdef ENABLE_SIM
	initcond_generate(mas, rad, pos, vel, N);
	#endif

	TIMING_STOP();
	TIMING_PRINT("simulation_init()", "initcond_generate");
	TIMING_START();

	#ifdef ENABLE_SIM
	solver_run(acc, pot, rho, prs, mas, rad, pos, vel, N, 0u);
	#endif

	TIMING_STOP();
	TIMING_PRINT("simulation_init()", "solver_run");
	TIMING_START();

	{
		sim.step_number = 0u;

		sim.N   = N;

		sim.mem = mem;

		sim.mas = mas;
		sim.rad = rad;
		sim.pos = pos;
		sim.vel = vel;

		sim.acc = acc;

		sim.pot = pot;

		sim.rho = rho;
		sim.prs = prs;

		sim.ken = ken;
		sim.pen = pen;
	}

	*simulation = sim;

	return STARFLOOD_SUCCESS;
}

int simulation_free(sim_t* restrict simulation) {
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

	sim_t sim = *simulation;

	TIMING_START();

	free(sim.mem);

	TIMING_STOP();
	TIMING_PRINT("simulation_free()", "free()");

	{
		sim.step_number = 0u;

		sim.N = 0u;

		sim.mem = NULL;

		sim.mas = (real*)NULL;
		sim.rad = (real*)NULL;
		sim.pos = (real*)NULL;
		sim.vel = (real*)NULL;

		sim.acc = (real*)NULL;

		sim.pot = (real*)NULL;

		sim.rho = (real*)NULL;
		sim.prs = (real*)NULL;

		sim.ken = (real*)NULL;
		sim.pen = (real*)NULL;
	}

	*simulation = sim;

	return STARFLOOD_SUCCESS;
}

int simulation_read(sim_t* restrict simulation, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *simulation;

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

	#ifdef ENABLE_SPH
	fread(sim.mem, sizeof(real), (size_t)N * (size_t)(1u+1u+3u+3u+3u+1u+1u+1u+1u+1u), file);
	#else
	fread(sim.mem, sizeof(real), (size_t)N * (size_t)(1u+3u+3u+3u+1u+1u+1u), file);
	#endif

	/*
	fread(sim.mas, sizeof(float), (size_t)N * (size_t)(1u), file);
	#ifdef ENABLE_SPH
	fread(sim.rad, sizeof(float), (size_t)N * (size_t)(1u), file);
	#endif
	fread(sim.pos, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.vel, sizeof(float), (size_t)N * (size_t)(3u), file);
	fread(sim.acc, sizeof(float), (size_t)N * (size_t)(3u), file);
	//fread(sim.pot, sizeof(float), (size_t)N * (size_t)(1u), file);
	//fread(sim.kin, sizeof(float), (size_t)N * (size_t)(1u), file);
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

int simulation_save(const sim_t* restrict simulation, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *simulation;

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

	#ifdef ENABLE_SPH
	fwrite(sim.mem, sizeof(real), (size_t)N * (size_t)(1u+1u+3u+3u+3u+1u+1u+1u+1u+1u), file);
	#else
	fwrite(sim.mem, sizeof(real), (size_t)N * (size_t)(1u+3u+3u+3u+1u+1u+1u), file);
	#endif

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

int simulation_step(sim_t* restrict simulation) {
	TIMING_INIT();

	sim_t sim = *simulation;

	const real dt = (real)TIMESTEP_SIZE;

	unsigned int step_number = sim.step_number;

	unsigned int N = sim.N;

	real* mas = sim.mas;
	real* rad = sim.rad;
	real* pos = sim.pos;
	real* vel = sim.vel;

	real* acc = sim.acc;

	real* pot = sim.pot;

	real* rho = sim.rho;
	real* prs = sim.prs;

	real* ken = sim.ken;
	real* pen = sim.pen;

	// for the very first leapfrog kick
	//if(0u >= step_number) {
	//	solver_run(acc, pot, mas, rad, pos, N, step_number);
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

	// run the solver
	solver_run(acc, pot, rho, prs, mas, rad, pos, vel, N, step_number);

	TIMING_STOP();
	TIMING_PRINT("simulation_step()", "solver_run");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif

	// update kinetic energy
	for(unsigned int i = 0u; i < N; i++) {
		real m_i = mas[i];

		real v_i[3] = {
			vel[3u*i+0u],
			vel[3u*i+1u],
			vel[3u*i+2u]
		};

		real v2 = (v_i[0]*v_i[0])+(v_i[1]*v_i[1])+(v_i[2]*v_i[2]);

		ken[i] = (real)0.5 * m_i * v2; // K = (1/2) * m * v^2
	}

	// update potential energy
	for(unsigned int i = 0u; i < N; i++) {
		pen[i] = mas[i] * pot[i];
	}

	// compute total energy
	{
		double T_sum = 0.0, T_c = 0.0; // Kinetic Energy
		double V_sum = 0.0, V_c = 0.0; // Potential Energy

		TIMING_START();

		for(unsigned int i = 0u; i < N; i++) {
			// Naïve summation
			//V_sum += (double)ken[i];

			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			double y = (double)ken[i] - T_c;
			volatile double t = T_sum + y;
			volatile double z = t - T_sum;
			T_c = z - y;
			T_sum = t;
		}

		TIMING_STOP();
		TIMING_PRINT("simulation_step()", "ken_sum");
		#ifdef LOG_TIMINGS_SIM_STEP
		LOG_TIMING(log_timings_sim_step);
		#endif
		TIMING_START();

		for(unsigned int i = 0u; i < N; i++) {
			// Naïve summation
			//V_sum += (double)pen[i];

			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			double y = (double)pen[i] - V_c;
			volatile double t = V_sum + y;
			volatile double z = t - V_sum;
			V_c = z - y;
			V_sum = t;
		}

		TIMING_STOP();
		TIMING_PRINT("simulation_step()", "pen_sum");
		#ifdef LOG_TIMINGS_SIM_STEP
		LOG_TIMING(log_timings_sim_step);
		#endif

		#ifdef LOG_STATISTICS
		fprintf(log_statistics.file, ",%.015f,%.015f,%.015f", V_sum + T_sum, V_sum, T_sum);
		#endif

		printf("E (tot) = % .015f\nV (pot) = % .015f\nT (kin) = % .015f\n", V_sum + T_sum, V_sum, T_sum);
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
