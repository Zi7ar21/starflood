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

// Maximum property name length
#define PLY_MAX_PROPERTY_LEN 32

// Maximum number of properties
#define PLY_MAX_PROPERTY_NUM 32

/*
// Enumeration of PLY format types
enum ply_format_t {
	PLY_FORMAT_ASCII,
	PLY_FORMAT_BINARY_LE,
	PLY_FORMAT_BINARY_BE
};
*/

// ply property type
typedef struct {
	char name[PLY_MAX_PROPERTY_LEN];
	enum sim_conf sim_parameter;
	real* data;
} ply_property_t;

// ply header type
typedef struct {
	//enum ply_format_t;
	size_t num_vertices;
	size_t num_property;
	ply_property_t property[PLY_MAX_PROPERTY_NUM];
} ply_header_t;

int sim_read(sim_t* restrict sim_ptr, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

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
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "sim_read()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "rb");

	TIMING_STOP();
	TIMING_PRINT("sim_read()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "sim_read()", filename, "rb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	#if ( (1 == SIM_FILE_FORMAT_READ) || (2 == SIM_FILE_FORMAT_READ) )
	if( 0 != fseek(file, 0l, SEEK_END) ) {
		fprintf(stderr, "%s error: fseek(file, 0l, SEEK_END) ", "sim_read()");
		perror("failed");
		fclose(file);
		return STARFLOOD_FAILURE;
	}

	long file_length = ftell(file);

	if(0l > file_length) {
		fprintf(stderr, "%s error: ftell(file) ", "sim_read()");
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
	TIMING_PRINT("sim_read()", "parse_header");
	TIMING_START();
	#endif

	size_t mem_length = (size_t)sim.N * (size_t)SIM_DOF;

	#if (0 == SIM_FILE_FORMAT_READ)
	
	if( fread(sim.mem, sizeof(real), mem_length, file) ) {
	}

	TIMING_STOP();
	TIMING_PRINT("sim_read()", "fread()");
	TIMING_START();
	#endif

	#if (1 == SIM_FILE_FORMAT_READ)
	fread(sim.pos, sizeof(real), (size_t)3u * (size_t)sim.N, file);

	TIMING_STOP();
	TIMING_PRINT("sim_read()", "fwrite()");
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
	TIMING_PRINT("sim_read()", "loop_fread");
	TIMING_START();
	#endif

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "sim_read()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_read()", "fclose()");

	return STARFLOOD_SUCCESS;
}

int sim_save(const sim_t* restrict sim_ptr, const char* restrict filename) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

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
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "sim_save()");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("sim_save()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "sim_save()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	#if ( (1 == SIM_FILE_FORMAT_SAVE) || (2 == SIM_FILE_FORMAT_SAVE) )
	ply_header_t ply_header;

	ply_header.num_vertices = (size_t)sim.N;

	#if (1 == SIM_FILE_FORMAT_SAVE)
	ply_header.num_property = (size_t)3u;
	#else
	ply_header.num_property = (size_t)SIM_DOF;
	#endif

	for(size_t i = (size_t)0u; i < ply_header.num_property; i++) {
		char property_name[PLY_MAX_PROPERTY_LEN];

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
		ply_header.property[i].sim_parameter = (enum sim_conf)i;
		ply_header.property[i].data = sim_find(&sim, (enum sim_conf)i);
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
	TIMING_PRINT("sim_save()", "print_header");
	TIMING_START();
	#endif

	if( 0 != fflush(file) ) {
		fprintf(stderr, "%s error: fflush(file) ", "sim_save()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save()", "fflush()");
	TIMING_START();

	#if (0 == SIM_FILE_FORMAT_SAVE)
	fwrite(sim.mem, sizeof(real), (size_t)sim.N * (size_t)SIM_DOF, file);

	TIMING_STOP();
	TIMING_PRINT("sim_save()", "fwrite()");
	TIMING_START();
	#endif

	#if (1 == SIM_FILE_FORMAT_SAVE)
	real* pos = sim_find(&sim, SIM_POS);

	fwrite(pos, sizeof(real), ply_header.num_vertices * ply_header.num_property, file);

	TIMING_STOP();
	TIMING_PRINT("sim_save()", "fwrite()");
	TIMING_START();
	#endif

	#if (2 == SIM_FILE_FORMAT_SAVE)
	for(unsigned int i = 0u; i < ply_header.num_vertices; i++) {
		for(unsigned int j = 0u; j < ply_header.num_property; j++) {
			enum sim_conf sim_parameter = ply_header.property[j].sim_parameter;

			real property = (real)0.0;

			if(9u > j) {
				property = ply_header.property[j].data[3u*i+(j%3u)];
			} else {
				property = ply_header.property[j].data[i];
			}

			fwrite(&property, sizeof(real), (size_t)1u, file);
		}
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save()", "loop_fwrite");
	TIMING_START();
	#endif

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "sim_save()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_save()", "fclose()");

	return STARFLOOD_SUCCESS;
}
