// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "simulation.h"

#include <stddef.h>
#include <stdio.h>

#include "common.h"
#include "config.h"
#include "timing.h"
#include "types.h"

struct ply_properties {
	// pos
	real x;
	real y;
	real z;

	real mas;

	#ifdef ENABLE_SPH
	real rad;
	#endif

	real vel_x;
	real vel_y;
	real vel_z;

	real acc_x;
	real acc_y;
	real acc_z;

	real pot;

	#ifdef ENABLE_SPH
	real rho;
	real prs;
	#endif

	real ken;
	real pen;
};

int simulation_read(sim_t* restrict simulation, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *simulation;

	#ifdef STARFLOOD_DOUBLE_PRECISION
	if( sizeof(real) != sizeof(double) ) {
		return STARFLOOD_FAILURE;
	}
	#else
	if( sizeof(real) != sizeof(float) ) {
		return STARFLOOD_FAILURE;
	}
	#endif

	if(NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	#ifdef STARFLOOD_DOUBLE_PRECISION
	volatile uint64_t magic_number = (uint32_t)0x401921FB54442D18u;

	uint64_t accumulator = (uint64_t)0u;

	for(size_t i = (size_t)0u; i < sizeof(magic_number); i++) {
		accumulator |= (uint64_t)( (uint8_t*)&magic_number )[i] << (uint64_t)(i * (size_t)8u);
	}

	int byte_order_is_le = (uint64_t)0x401921FB54442D18u == accumulator;
	int byte_order_is_be = (uint64_t)0x182D4454FB211940u == accumulator;
	#else
	volatile uint32_t magic_number = (uint32_t)0x40C90FDBu;

	uint32_t accumulator = (uint32_t)0u;

	for(size_t i = (size_t)0u; i < sizeof(magic_number); i++) {
		accumulator |= (uint32_t)( (uint8_t*)&magic_number )[i] << (uint32_t)(i * (size_t)8u);
	}

	int byte_order_is_le = (uint32_t)0x40C90FDBu == accumulator;
	int byte_order_is_be = (uint32_t)0xDB0FC940u == accumulator;
	#endif

	if( (!byte_order_is_le && !byte_order_is_be) || (byte_order_is_be && byte_order_is_le) ) {
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "simulation_read()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "rb");

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "simulation_read()", filename, "rb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	#if ( (1 == SIM_FILE_FORMAT_READ) || (2 == SIM_FILE_FORMAT_READ) )
	if( 0 != fseek(file, 0l, SEEK_END) ) {
		fprintf(stderr, "%s error: fseek(file, 0l, SEEK_END) ", "simulation_read()");
		perror("failed");
		fclose(file);
		return STARFLOOD_FAILURE;
	}

	long file_length = ftell(file);

	if(0l > file_length) {
		fprintf(stderr, "%s error: ftell(file) ", "simulation_read()");
		perror("failed");
		fclose(file);
		return STARFLOOD_FAILURE;
	}

	//fread();

	/*
	fprintf(file,
		"ply\n"
		"format %s 1.0\n"
		"comment Created by Starflood version %d.%d.%d\n"
		"element vertex %zu\n"
		"property " REAL_TYPE_STRING " x\n"
		"property " REAL_TYPE_STRING " y\n"
		"property " REAL_TYPE_STRING " z\n"
		#if (2 == SIM_FILE_FORMAT_SAVE)
		"property " REAL_TYPE_STRING " mas\n"
		#ifdef ENABLE_SPH
		"property " REAL_TYPE_STRING " rad\n"
		#endif
		"property " REAL_TYPE_STRING " vel_x\n"
		"property " REAL_TYPE_STRING " vel_y\n"
		"property " REAL_TYPE_STRING " vel_z\n"
		"property " REAL_TYPE_STRING " acc_x\n"
		"property " REAL_TYPE_STRING " acc_y\n"
		"property " REAL_TYPE_STRING " acc_z\n"
		"property " REAL_TYPE_STRING " pot\n"
		#ifdef ENABLE_SPH
		"property " REAL_TYPE_STRING " rho\n"
		"property " REAL_TYPE_STRING " prs\n"
		#endif
		"property " REAL_TYPE_STRING " pen\n"
		"property " REAL_TYPE_STRING " ken\n"
		#endif
		"end_header\n",
		(byte_order_is_le ? "binary_little_endian" : "binary_big_endian"),
		STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH,
		(size_t)sim.N
	);
	*/

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "parse_header");
	TIMING_START();
	#endif

	#ifdef ENABLE_SPH
	size_t mem_length = (size_t)sim.N * (size_t)(1u+0u+3u+3u+3u+1u+0u+0u+1u+1u);
	#else
	size_t mem_length = (size_t)sim.N * (size_t)(1u+1u+3u+3u+3u+1u+1u+1u+1u+1u);
	#endif

	#if (0 == SIM_FILE_FORMAT_READ)
	
	if( fread(sim.mem, sizeof(real), mem_length, file) ) {
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fread()");
	TIMING_START();
	#endif

	#if (1 == SIM_FILE_FORMAT_READ)
	fread(sim.pos, sizeof(real), (size_t)3u * (size_t)sim.N, file);

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fwrite()");
	TIMING_START();
	#endif

	#if (2 == SIM_FILE_FORMAT_READ)
	for(unsigned int i = 0u; i < sim.N; i++) {
		fread(&sim.pos[3u*i], sizeof(real), (size_t)3u, file);
		fread(&sim.mas[   i], sizeof(real), (size_t)1u, file);
		#ifdef ENABLE_SPH
		fread(&sim.rad[   i], sizeof(real), (size_t)1u, file);
		#endif
		fread(&sim.vel[3u*i], sizeof(real), (size_t)3u, file);
		fread(&sim.acc[3u*i], sizeof(real), (size_t)3u, file);
		fread(&sim.pot[   i], sizeof(real), (size_t)1u, file);
		#ifdef ENABLE_SPH
		fread(&sim.rho[   i], sizeof(real), (size_t)1u, file);
		fread(&sim.prs[   i], sizeof(real), (size_t)1u, file);
		#endif
		fread(&sim.ken[   i], sizeof(real), (size_t)1u, file);
		fread(&sim.pen[   i], sizeof(real), (size_t)1u, file);
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "loop_fread");
	TIMING_START();
	#endif

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "simulation_read()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_read()", "fclose()");

	return STARFLOOD_SUCCESS;
}

int simulation_save(const sim_t* restrict simulation, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *simulation;

	#ifdef STARFLOOD_DOUBLE_PRECISION
	if( sizeof(real) != sizeof(double) ) {
		return STARFLOOD_FAILURE;
	}
	#else
	if( sizeof(real) != sizeof(float) ) {
		return STARFLOOD_FAILURE;
	}
	#endif

	if(NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	#ifdef STARFLOOD_DOUBLE_PRECISION
	volatile uint64_t magic_number = (uint32_t)0x401921FB54442D18u;

	uint64_t accumulator = (uint64_t)0u;

	for(size_t i = (size_t)0u; i < sizeof(magic_number); i++) {
		accumulator |= (uint64_t)( (uint8_t*)&magic_number )[i] << (uint64_t)(i * (size_t)8u);
	}

	int byte_order_is_le = (uint64_t)0x401921FB54442D18u == accumulator;
	int byte_order_is_be = (uint64_t)0x182D4454FB211940u == accumulator;
	#else
	volatile uint32_t magic_number = (uint32_t)0x40C90FDBu;

	uint32_t accumulator = (uint32_t)0u;

	for(size_t i = (size_t)0u; i < sizeof(magic_number); i++) {
		accumulator |= (uint32_t)( (uint8_t*)&magic_number )[i] << (uint32_t)(i * (size_t)8u);
	}

	int byte_order_is_le = (uint32_t)0x40C90FDBu == accumulator;
	int byte_order_is_be = (uint32_t)0xDB0FC940u == accumulator;
	#endif

	if( (!byte_order_is_le && !byte_order_is_be) || (byte_order_is_be && byte_order_is_le) ) {
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "simulation_save()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "simulation_save()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	#if ( (1 == SIM_FILE_FORMAT_SAVE) || (2 == SIM_FILE_FORMAT_SAVE) )
	fprintf(file,
		"ply\n"
		"format %s 1.0\n"
		"comment Created by Starflood version %d.%d.%d\n"
		"element vertex %zu\n"
		"property " REAL_TYPE_STRING " x\n"
		"property " REAL_TYPE_STRING " y\n"
		"property " REAL_TYPE_STRING " z\n"
		#if (2 == SIM_FILE_FORMAT_SAVE)
		"property " REAL_TYPE_STRING " mas\n"
		#ifdef ENABLE_SPH
		"property " REAL_TYPE_STRING " rad\n"
		#endif
		"property " REAL_TYPE_STRING " vel_x\n"
		"property " REAL_TYPE_STRING " vel_y\n"
		"property " REAL_TYPE_STRING " vel_z\n"
		"property " REAL_TYPE_STRING " acc_x\n"
		"property " REAL_TYPE_STRING " acc_y\n"
		"property " REAL_TYPE_STRING " acc_z\n"
		"property " REAL_TYPE_STRING " pot\n"
		#ifdef ENABLE_SPH
		"property " REAL_TYPE_STRING " rho\n"
		"property " REAL_TYPE_STRING " prs\n"
		#endif
		"property " REAL_TYPE_STRING " pen\n"
		"property " REAL_TYPE_STRING " ken\n"
		#endif
		"end_header\n",
		(byte_order_is_le ? "binary_little_endian" : "binary_big_endian"),
		STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH,
		(size_t)sim.N
	);

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "print_header");
	TIMING_START();
	#endif

	if( 0 != fflush(file) ) {
		fprintf(stderr, "%s error: fflush(file) ", "simulation_save()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fflush()");
	TIMING_START();

	#if (0 == SIM_FILE_FORMAT_SAVE)
	#ifdef ENABLE_SPH
	fwrite(sim.mem, sizeof(real), (size_t)sim.N * (size_t)(1u+0u+3u+3u+3u+1u+0u+0u+1u+1u), file);
	#else
	fwrite(sim.mem, sizeof(real), (size_t)sim.N * (size_t)(1u+1u+3u+3u+3u+1u+1u+1u+1u+1u), file);
	#endif

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fwrite()");
	TIMING_START();
	#endif

	#if (1 == SIM_FILE_FORMAT_SAVE)
	fwrite(sim.pos, sizeof(real), (size_t)3u * (size_t)sim.N, file);

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fwrite()");
	TIMING_START();
	#endif

	#if (2 == SIM_FILE_FORMAT_SAVE)
	for(unsigned int i = 0u; i < sim.N; i++) {
		struct ply_properties ply_prop = {
			sim.pos[3u*i+0u],
			sim.pos[3u*i+1u],
			sim.pos[3u*i+2u],
			sim.mas[   i   ],
			#ifdef ENABLE_SPH
			sim.rad[   i   ],
			#endif
			sim.vel[3u*i+0u],
			sim.vel[3u*i+1u],
			sim.vel[3u*i+2u],
			sim.acc[3u*i+0u],
			sim.acc[3u*i+1u],
			sim.acc[3u*i+2u],
			sim.pot[   i   ],
			#ifdef ENABLE_SPH
			sim.rho[   i   ],
			sim.prs[   i   ],
			#endif
			sim.ken[   i   ],
			sim.pen[   i   ]
		};

		fwrite(&ply_prop, sizeof(ply_prop), (size_t)1u, file);
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "loop_fwrite");
	TIMING_START();
	#endif

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "simulation_save()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("simulation_save()", "fclose()");

	return STARFLOOD_SUCCESS;
}
