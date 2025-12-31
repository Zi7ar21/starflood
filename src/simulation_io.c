// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "simulation.h"

#include <stddef.h>
#include <stdio.h>

#include "common.h"
#include "config.h"
#include "timing.h"
#include "types.h"

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
