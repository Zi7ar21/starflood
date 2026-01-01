#include "ply.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "common.h"
#include "types.h"

int ply_write(const char* restrict filename, real* restrict vertex_xyz, size_t num_vertices) {
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

	FILE* ply_file = fopen(filename, "wb");

	if(NULL == (void*)ply_file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "ply_write()", filename, "wb");
		perror("failed");
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
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "ply_write()");
		return STARFLOOD_FAILURE;
	}

	fprintf(ply_file,
		"ply\n"
		"format %s 1.0\n"
		"comment Created by Starflood version %d.%d.%d\n"
		"element vertex %zu\n"
		#ifdef STARFLOOD_DOUBLE_PRECISION
		"property float x\n"
		"property float y\n"
		"property float z\n"
		#else
		"property double x\n"
		"property double y\n"
		"property double z\n"
		#endif
		"end_header\n",
		(byte_order_is_le ? "binary_little_endian" : "binary_big_endian"),
		STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH,
		num_vertices
	);

	if( 0 != fflush(ply_file) ) {
		fprintf(stderr, "%s error: fflush() ", "ply_write()");
		perror("failed");
	}

	fwrite(vertex_xyz, sizeof(float), (size_t)3u * num_vertices, ply_file);

	if( 0 != fclose(ply_file) ) {
		fprintf(stderr, "%s error: fclose() ", "ply_write()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}
