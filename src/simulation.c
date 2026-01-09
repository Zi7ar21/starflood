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
#include "grid.h"
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

// Locates a simulation parameter
real* sim_find(const sim_t* restrict sim_ptr, enum sim_conf parameter_enum) {
	sim_t sim = *sim_ptr;

	if(SIM_DOF <= parameter_enum) {
		return (real*)NULL;
	}

	if(_SIM_POS_0 == parameter_enum || _SIM_POS_1 == parameter_enum || _SIM_POS_2 == parameter_enum) {
		parameter_enum = SIM_POS;
	}

	if(_SIM_VEL_0 == parameter_enum || _SIM_VEL_1 == parameter_enum || _SIM_VEL_2 == parameter_enum) {
		parameter_enum = SIM_VEL;
	}

	if(_SIM_ACC_0 == parameter_enum || _SIM_ACC_1 == parameter_enum || _SIM_ACC_2 == parameter_enum) {
		parameter_enum = SIM_ACC;
	}

	return &( (real*)sim.mem )[(size_t)sim.N * (size_t)parameter_enum];
}

int sim_init(sim_t* restrict sim_ptr, unsigned int N) {
	TIMING_INIT();

	#ifdef LOG_STATISTICS
	if( STARFLOOD_SUCCESS != log_init(&log_statistics, OUTPUT_DIR "/" LOG_STATISTICS) ) {
		fprintf(stderr, "%s error: %s failed!\n", "sim_init()", "log_init(&log_statistics)");
		return STARFLOOD_FAILURE;
	}

	fprintf(log_statistics.file, "%s,%s,%s,%s\n", "step_number", "total_energy", "potential_energy", "kinetic_energy");
	log_sync(&log_statistics);
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	if( STARFLOOD_SUCCESS != log_init(&log_timings_sim_step, OUTPUT_DIR "/" LOG_TIMINGS_SIM_STEP) ) {
		fprintf(stderr, "%s error: %s failed!\n", "sim_init()", "log_init(&log_timings_sim_step)");
		return STARFLOOD_FAILURE;
	}

	fprintf(log_timings_sim_step.file, "%s,%s,%s,%s,%s,%s,%s\n", "step_number", "kick_0", "drift", "solver_run", "kick_1", "ken_sum", "pen_sum");
	log_sync(&log_timings_sim_step);
	#endif

	sim_t sim = *sim_ptr;

	void* mem = NULL;

	size_t mem_size = sizeof(real) * (size_t)SIM_DOF * (size_t)N;

	printf("Simulation Memory Addresses:\n");

	TIMING_START();

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
	#endif

	TIMING_STOP();

	#ifdef STARFLOOD_ALIGNMENT
	TIMING_PRINT("sim_init()", "posix_memalign()");
	#else
	TIMING_PRINT("sim_init()", "malloc()");
	#endif

	if(NULL == mem) {
		#ifdef STARFLOOD_ALIGNMENT
		fprintf(stderr, "%s error: mem is NULL after posix_memalign(&mem, %zu, %zu", "sim_init()", (size_t)STARFLOOD_ALIGNMENT, mem_size);
		#else
		fprintf(stderr, "%s error: mem is NULL after malloc(%zu", "sim_init()", mem_size);
		#endif

		perror(")");

		return STARFLOOD_FAILURE;
	}

	printf("  mem: %p (%zu bytes)\n", mem, mem_size);

	TIMING_START();

	memset(mem, 0, mem_size);

	TIMING_STOP();
	TIMING_PRINT("sim_init()", "memset()");

	/*
	for(int i = 0; i < SIM_DOF; i++) {
		if(i 
		switch(i) {
			case(SIM_POS) {
			}
		}
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
	*/

	#ifdef ENABLE_GRID
	TIMING_START();

	real grid_bounds_min[3] = {
		(real)(-2.000e0),
		(real)(-2.000e0),
		(real)(-2.000e0)
	};

	real grid_bounds_max[3] = {
		(real)( 2.000e0),
		(real)( 2.000e0),
		(real)( 2.000e0)
	};

	if( STARFLOOD_SUCCESS != grid_init(&sim.grid, 1u << (unsigned int)GRID_SIZE, grid_bounds_min, grid_bounds_max) ) {
	}

	TIMING_STOP();
	TIMING_PRINT("sim_init()", "grid_init()");
	#endif

	{
		sim.step_number = 0u;

		sim.N = N;

		sim.mem = mem;
	}

	*sim_ptr = sim;

	return STARFLOOD_SUCCESS;
}

int sim_free(sim_t* restrict sim_ptr) {
	TIMING_INIT();

	#ifdef LOG_STATISTICS
	if( STARFLOOD_SUCCESS != log_free(&log_statistics) ) {
		fprintf(stderr, "%s error: %s failed!\n", "sim_free()", "log_free(&log_statistics)");
	}
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	if( STARFLOOD_SUCCESS != log_free(&log_timings_sim_step) ) {
		fprintf(stderr, "%s error: %s failed!\n", "sim_free()", "log_free(&log_timings_sim_step)");
	}
	#endif

	sim_t sim = *sim_ptr;

	TIMING_START();

	free(sim.mem);

	TIMING_STOP();
	TIMING_PRINT("sim_free()", "free()");

	#ifdef ENABLE_GRID
	grid_free(&sim.grid);
	#endif

	{
		sim.step_number = 0u;

		sim.N = 0u;

		sim.mem = NULL;
	}

	*sim_ptr = sim;

	return STARFLOOD_SUCCESS;
}

int sim_step(sim_t* restrict sim_ptr) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	unsigned int step_number = sim.step_number;
	unsigned int N = sim.N;
	real* pos = sim_find(&sim, SIM_POS);
	real* vel = sim_find(&sim, SIM_VEL);
	real* acc = sim_find(&sim, SIM_ACC);
	real* mas = sim_find(&sim, SIM_MAS);
	real* pot = sim_find(&sim, SIM_POT);
	real* ken = sim_find(&sim, SIM_KEN);
	real* pen = sim_find(&sim, SIM_PEN);

	const real dt = (real)TIMESTEP_SIZE;

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
	TIMING_PRINT("sim_step()", "kick_0");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif
	TIMING_START();

	// drift
	for(unsigned int i = 0u; i < 3u * N; i++) {
		pos[i] += dt * vel[i];
	}

	TIMING_STOP();
	TIMING_PRINT("sim_step()", "drift");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif
	TIMING_START();

	// run the solver
	if( STARFLOOD_SUCCESS != solver_run(&sim, step_number) ) {
		fprintf(stderr, "%s error: %s failed.\n", "sim_step()", "solver_run()");
	}

	TIMING_STOP();
	TIMING_PRINT("sim_step()", "solver_run");
	#ifdef LOG_TIMINGS_SIM_STEP
	LOG_TIMING(log_timings_sim_step);
	#endif
	TIMING_START();

	// kick_1
	for(unsigned int i = 0u; i < 3u * N; i++) {
		vel[i] += (real)0.5 * dt * acc[i];
	}

	TIMING_STOP();	
	TIMING_PRINT("sim_step()", "kick_1");
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
		double T_sum = 0.0, T_com = 0.0; // Kinetic Energy
		double V_sum = 0.0, V_com = 0.0; // Potential Energy

		TIMING_START();

		for(unsigned int i = 0u; i < N; i++) {
			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			double y = (double)ken[i] - T_com;
			volatile double t = T_sum + y;
			volatile double z = t - T_sum;
			T_com = z - y;
			T_sum = t;

			// Naïve summation
			//T_sum += (double)ken[i];
		}

		TIMING_STOP();
		TIMING_PRINT("sim_step()", "ken_sum");
		#ifdef LOG_TIMINGS_SIM_STEP
		LOG_TIMING(log_timings_sim_step);
		#endif
		TIMING_START();

		for(unsigned int i = 0u; i < N; i++) {
			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			double y = (double)pen[i] - V_com;
			volatile double t = V_sum + y;
			volatile double z = t - V_sum;
			V_com = z - y;
			V_sum = t;

			// Naïve summation
			//V_sum += (double)pen[i];
		}

		TIMING_STOP();
		TIMING_PRINT("sim_step()", "pen_sum");
		#ifdef LOG_TIMINGS_SIM_STEP
		LOG_TIMING(log_timings_sim_step);
		#endif

		#ifdef LOG_STATISTICS
		fprintf(log_statistics.file, ",%.015f,%.015f,%.015f", V_sum + T_sum, V_sum, T_sum);
		#endif

		printf("E (tot) = % .015f\nV (pot) = % .015f\nT (kin) = % .015f\n", V_sum + T_sum, V_sum, T_sum);
	}

	#ifdef LOG_STATISTICS
	fprintf(log_statistics.file, "\n");
	log_sync(&log_statistics);
	#endif

	#ifdef LOG_TIMINGS_SIM_STEP
	fprintf(log_timings_sim_step.file, "\n");
	log_sync(&log_timings_sim_step);
	#endif

	*sim_ptr = sim;

	return STARFLOOD_SUCCESS;
}
