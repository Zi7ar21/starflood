// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "simulation.h"

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "common.h"
#include "config.h"
#include "timing.h"
#include "types.h"

int sim_read_raw(sim_t* restrict sim_ptr, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	/*
	#ifdef STARFLOOD_DOUBLE_PRECISION
	if( sizeof(real) != sizeof(f64) ) {
		return STARFLOOD_FAILURE;
	}
	#else
	if( sizeof(real) != sizeof(f32) ) {
		return STARFLOOD_FAILURE;
	}
	#endif
	*/

	if(NULL == (void*)filename) {
		fprintf(stderr, "%s error: filename is NULL!\n", "sim_read_raw()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "rb");

	TIMING_STOP();
	TIMING_PRINT("sim_read_raw()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "sim_read_raw()", filename, "rb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	TIMING_STOP();
	TIMING_PRINT("sim_read_raw()", "parse_header");

	const size_t mem_length = (size_t)sim.N * (size_t)SIM_DOF;

	TIMING_START();

	if( mem_length != fread(sim.mem, sizeof(real), mem_length, file) ) {
		fprintf(stderr, "%s error: fread(file) ", "sim_read_raw()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("sim_read_raw()", "fread()");
	TIMING_START();

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "sim_read_raw()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_read_raw()", "fclose()");

	return STARFLOOD_SUCCESS;
}

int sim_save_raw(const sim_t* restrict sim_ptr, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	/*
	#ifdef STARFLOOD_DOUBLE_PRECISION
	if( sizeof(real) != sizeof(f64) ) {
		return STARFLOOD_FAILURE;
	}
	#else
	if( sizeof(real) != sizeof(f32) ) {
		return STARFLOOD_FAILURE;
	}
	#endif
	*/

	if(NULL == (void*)filename) {
		fprintf(stderr, "%s error: filename is NULL!\n", "sim_save_raw()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("sim_save_raw()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "sim_save_raw()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	if( 0 != fflush(file) ) {
		fprintf(stderr, "%s error: fflush(file) ", "sim_save_raw()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save_raw()", "fflush()");

	const size_t mem_length = (size_t)sim.N * (size_t)SIM_DOF;

	TIMING_START();

	if( mem_length != fwrite(sim.mem, sizeof(real), mem_length, file) ) {
		fprintf(stderr, "%s error: fwrite(file) ", "sim_save_raw()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save_raw()", "fwrite()");
	TIMING_START();

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "sim_save_raw()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save_raw()", "fclose()");

	return STARFLOOD_SUCCESS;
}

// Maximum property name length
#define IO_PLY_MAX_PROPERTY_LEN 32

// Maximum number of properties
#define IO_PLY_MAX_PROPERTY_NUM 32

// Enumeration of format types
enum io_ply_format_t {
	IO_PLY_FORMAT_ASCII,
	IO_PLY_FORMAT_BINARY_LE,
	IO_PLY_FORMAT_BINARY_BE
};

typedef struct {
	enum sim_conf sim_parameter;
	char name[IO_PLY_MAX_PROPERTY_LEN];
} io_ply_property_t;

typedef struct {
	enum io_ply_format_t format;
	size_t num_vertices;
	size_t num_property;
	io_ply_property_t property[IO_PLY_MAX_PROPERTY_NUM];
} io_ply_header_t;

int sim_read_ply(sim_t* restrict sim_ptr, const char* restrict filename) {
	TIMING_INIT();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	if( sizeof(real) != sizeof(f64) ) {
		return STARFLOOD_FAILURE;
	}
	#else
	if( sizeof(real) != sizeof(f32) ) {
		return STARFLOOD_FAILURE;
	}
	#endif

	if(NULL == (void*)filename) {
		fprintf(stderr, "%s error: filename is NULL!\n", "sim_read_ply()");
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

	return STARFLOOD_FAILURE;
}

int sim_save_ply(const sim_t* restrict sim_ptr, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	#ifdef STARFLOOD_DOUBLE_PRECISION
	if( sizeof(real) != sizeof(f64) ) {
		return STARFLOOD_FAILURE;
	}
	#else
	if( sizeof(real) != sizeof(f32) ) {
		return STARFLOOD_FAILURE;
	}
	#endif

	if(NULL == (void*)filename) {
		fprintf(stderr, "%s error: filename is NULL!\n", "sim_save_ply()");
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
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "sim_save_ply()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("sim_save_ply()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "sim_save_ply()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	io_ply_header_t ply_header;

	ply_header.num_vertices = (size_t)sim.N;

	#if (1 == SIM_FILE_FORMAT_SAVE)
	ply_header.num_property = (size_t)3u;
	#else
	ply_header.num_property = (size_t)SIM_DOF;
	#endif

	for(size_t i = (size_t)0u; i < ply_header.num_property; i++) {
		ply_header.property[i].sim_parameter = (enum sim_conf)i;

		char property_name[IO_PLY_MAX_PROPERTY_LEN];

		switch(i) {
			case _SIM_POS_0:
				strcpy(property_name, "x");
				break;
			case _SIM_POS_1:
				strcpy(property_name, "y");
				break;
			case _SIM_POS_2:
				strcpy(property_name, "z");
				break;
			case _SIM_VEL_0:
				strcpy(property_name, "vel_x");
				break;
			case _SIM_VEL_1:
				strcpy(property_name, "vel_y");
				break;
			case _SIM_VEL_2:
				strcpy(property_name, "vel_z");
				break;
			case _SIM_ACC_0:
				strcpy(property_name, "acc_x");
				break;
			case _SIM_ACC_1:
				strcpy(property_name, "acc_y");
				break;
			case _SIM_ACC_2:
				strcpy(property_name, "acc_z");
				break;
			case SIM_MAS:
				strcpy(property_name, "mas");
				break;
			case SIM_POT:
				strcpy(property_name, "pot");
				break;
			case SIM_KEN:
				strcpy(property_name, "ken");
				break;
			case SIM_PEN:
				strcpy(property_name, "pen");
				break;
			#ifdef ENABLE_SPH
			case SIM_RAD:
				strcpy(property_name, "rad");
				break;
			case SIM_RHO:
				strcpy(property_name, "rho");
				break;
			case SIM_PRS:
				strcpy(property_name, "prs");
				break;
			#endif
			default:
				break;
		}

		strcpy(ply_header.property[i].name, property_name);
	}

	fprintf(file,
		"ply\n"
		"format %s 1.0\n"
		"comment Created by Starflood version %d.%d.%d\n"
		"element vertex %zu\n",
		(byte_order_is_le ? "binary_little_endian" : "binary_big_endian"),
		STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH,
		(size_t)sim.N
	);

	for(size_t i = (size_t)0u; i < ply_header.num_property; i++) {
		fprintf(file, "property " REAL_TYPE_STRING " %s\n", ply_header.property[i].name);
	}

	fprintf(file, "end_header\n");

	TIMING_STOP();
	TIMING_PRINT("sim_save_ply()", "print_header");
	TIMING_START();

	if( 0 != fflush(file) ) {
		fprintf(stderr, "%s error: fflush(file) ", "sim_save_ply()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save_ply()", "fflush()");
	TIMING_START();

	#if (1 == SIM_FILE_FORMAT_SAVE)
	real* pos = sim_find(&sim, SIM_POS);

	fwrite(pos, sizeof(real), (size_t)3u * ply_header.num_vertices, file);

	TIMING_STOP();
	TIMING_PRINT("sim_save_ply()", "fwrite()");
	TIMING_START();
	#else
	for(unsigned int i = 0u; i < ply_header.num_vertices; i++) {
		for(unsigned int j = 0u; j < ply_header.num_property; j++) {
			enum sim_conf sim_parameter = ply_header.property[j].sim_parameter;

			real* data = sim_find(&sim, sim_parameter);

			real property = 9u > j ? data[3u * i + (j % 3u)] : data[i];

			fwrite(&property, sizeof(real), (size_t)1u, file);
		}
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save_ply()", "loop_fwrite");
	TIMING_START();
	#endif

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "sim_save_ply()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save_ply()", "fclose()");

	return STARFLOOD_SUCCESS;
}
