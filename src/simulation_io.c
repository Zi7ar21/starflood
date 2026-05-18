// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "simulation.h"

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "byte_order.h"
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

#define IO_PLY_HEADER_MAX_LINE_LEN 256
#define IO_PLY_HEADER_MAX_LINE_NUM 256

typedef enum {
	IO_PLY_HEADER_LINE_PLY,
	IO_PLY_HEADER_LINE_COMMENT,
	IO_PLY_HEADER_LINE_FORMAT,
	IO_PLY_HEADER_LINE_ELEMENT,
	IO_PLY_HEADER_LINE_PROPERTY,
	IO_PLY_HEADER_LINE_END,
} io_ply_header_line_t;

// Enumeration of format types
typedef enum {
	IO_PLY_FORMAT_ASCII,
	IO_PLY_FORMAT_BINARY_LE,
	IO_PLY_FORMAT_BINARY_BE
} io_ply_format_t;

typedef enum {
	IO_PLY_DATATYPE_FLOAT,
	IO_PLY_DATATYPE_DOUBLE,
} io_ply_datatype_t;

typedef struct {
	io_ply_datatype_t datatype;
	char name[IO_PLY_MAX_PROPERTY_LEN];
	enum sim_conf sim_parameter;
} io_ply_property_t;

typedef struct {
	io_ply_format_t format;
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

	fprintf(stderr, "%s error: PLY file reader implementation is incomplete.\n", "sim_read_ply()");

	/*
	if(NULL == (void*)filename) {
		fprintf(stderr, "%s error: filename is NULL!\n", "sim_read_ply()");
		return STARFLOOD_FAILURE;
	}

	byte_order_t host_byte_order = detect_host_byte_order();

	if(BYTE_ORDER_LE != host_byte_order && BYTE_ORDER_BE != host_byte_order) {
		fprintf(stderr, "%s error: PLY file I/O only works on machines standard byte ordering (little/big endian), but detect_host_byte_order returned %d.\n", "sim_read_ply()", (int)host_byte_order);
		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("sim_read_ply()", "fopen()");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "sim_read_ply()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	if( 0 != fseek(file, 0l, SEEK_END) ) {
		fprintf(stderr, "%s error: fseek(file, 0l, SEEK_END) ", "sim_read_ply()");
		perror("failed");
		fclose(file);
		return STARFLOOD_FAILURE;
	}

	long file_size = ftell(file);

	if(0l >= file_size) {
	}

	rewind(file);

	io_ply_header_t ply_header;

	int num_ply = 0;
	int num_form = 0;
	int num_elem = 0;
	int num_end = 0;

	TIMING_START();

	for(int line_num = 0; line_num < IO_PLY_HEADER_MAX_LINE_NUM; line_num++) {
		char header_line[IO_PLY_HEADER_MAX_LINE_LEN];

		fgets(header_line, IO_PLY_HEADER_MAX_LINE_LEN, file);

		if( ferror(file) ) {
			fprintf(stderr, "%s error: %s(file) set the file stream %s ", "sim_read_ply()", "ferror", "error");
			perror("indicator");
			fclose(file);
			return STARFLOOD_FAILURE;
		}

		if( feof(file) ) {
			fprintf(stderr, "%s error: %s(file) set the file stream %s ", "sim_read_ply()", "feof", "end");
			perror("indicator");
			fclose(file);
			return STARFLOOD_FAILURE;
		}

		char header_line_keyword[11];

		// "end_header" () is the longest first word of a line
		if( 1 != sscanf(header_line, "%10s", header_line_keyword) ) {
			fprintf(stderr, "%s error: sscanf(..., header_line_keyword) ", "sim_read_ply()");
			perror("failed");
			fclose(file);
			return STARFLOOD_FAILURE;
		}

		if( 0 == strcmp(header_line, "ply") ) {
			num_ply++;

			if(1 < num_ply) {
				fprintf(stderr, "%s error parsing PLY file header: \"ply\" found more than once.\n", "sim_read_ply()");
				fclose(file);
				return STARFLOOD_FAILURE;
			}

			continue;
		}

		if(1 > num_ply) {
			fprintf(stderr, "%s error parsing PLY file header: \"ply\" not found at the beginning of the file.\n", "sim_read_ply()");
			fclose(file);
			return STARFLOOD_FAILURE;
		}

		if( 0 == strcmp(header_line, "comment") ) {
			continue;
		}

		if( 0 == strcmp(header_line, "format") ) {
			num_form++;

			if(1 < num_form) {
				fprintf(stderr, "%s error parsing PLY file header: \"format\" found more than once.\n", "sim_read_ply()");
				fclose(file);
				return STARFLOOD_FAILURE;
			}

			continue;
		}

		if( 0 == strcmp(header_line, "element") ) {
			num_elem++;

			if(1 < num_elem) {
				fprintf(stderr, "%s error reading PLY file: Support for PLY files with multiple elements is not implemented yet.\n", "sim_read_ply()");
				fclose(file);
				return STARFLOOD_FAILURE;
			}

			char element_type[7];

			if( 1 != sscanf(header_line, "%*10s %6s", element_type) ) {
				fclose(file);
				return STARFLOOD_FAILURE;
			}

			unsigned int vertex_count;

			if( 0 != strcmp(element_type, "vertex") ) {
				fprintf(stderr, "%s error reading PLY file: Support for PLY files with non-vertex elements implemented yet.\n", "sim_read_ply()");
				fclose(file);
				return STARFLOOD_FAILURE;
			}

			continue;
		}

		if( 0 == strcmp(header_line, "property") ) {
			continue;
		}

		if( 0 == strcmp(header_line, "end_header") ) {
			num_end++;
			break;
		}

		fprintf(stderr, "%s error parsing PLY file header: Error while parsing line #%ld.\n", "sim_read_ply()", line_num);
		fclose(file);
		return STARFLOOD_FAILURE;
	}

	if(1 != num_end) {
		fprintf(stderr, "%s error parsing PLY file header: Exceeded max number of header lines (%d).\n", "sim_read_ply()", IO_PLY_HEADER_MAX_LINE_NUM);
		fclose(file);
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_read_ply()", "parse_header");
	TIMING_START();

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose(file) ", "sim_read_ply()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("sim_read_ply()", "fclose()");
	*/

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

	byte_order_t host_byte_order = detect_host_byte_order();

	if(BYTE_ORDER_LE != host_byte_order && BYTE_ORDER_BE != host_byte_order) {
		fprintf(stderr, "%s error: PLY file I/O only works on machines standard byte ordering (little/big endian), but detect_host_byte_order returned %d.\n", "sim_save_ply()", (int)host_byte_order);
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
		(BYTE_ORDER_LE == host_byte_order ? "binary_little_endian" : "binary_big_endian"),
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
