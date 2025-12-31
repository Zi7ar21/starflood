#include "ply.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "common.h"
#include "types.h"

int ply_write(const char* restrict filename, float* restrict vertex_xyz, size_t num_vertices) {
	if(NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	FILE* ply_file = fopen(filename, "wb");

	if(NULL == (void*)ply_file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "ply_write()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	volatile uint32_t magic_number = (uint32_t)0x40C90FDBu;
	//volatile uint64_t magic_number = (uint32_t)0x401921FB54442D18u;
	
	int byte_order_is_be = 1;
	int byte_order_is_le = 1;

	// This may seem a bit verbose, but it can handle:
	// - Arbitrary integer type sizes up to 2048-bit (replace uint32_t or uint64_t, and magic_number)
	// - Detection of non-little and non-big endianess (such as mixed endianess)
	// - Different endianess at runtime
	// - Machines with non-octet unsigned char (assuming uint8_t syntax works)
	for(size_t i = 0u; i < sizeof(magic_number); i++) {
		byte_order_is_le = ((uint8_t*)&magic_number)[i] == (uint8_t)((uint32_t)0xFFu & (magic_number >> (uint32_t)((                                  i)*(size_t)8u))) ? byte_order_is_le : 0;
		byte_order_is_be = ((uint8_t*)&magic_number)[i] == (uint8_t)((uint32_t)0xFFu & (magic_number >> (uint32_t)(((sizeof(magic_number)-(size_t)1u)-i)*(size_t)8u))) ? byte_order_is_be : 0;
		//byte_order_is_le = ((uint8_t*)&magic_number)[i] == (uint8_t)((uint64_t)0xFFu & (magic_number >> (uint64_t)((                                  i)*(size_t)8u))) ? byte_order_is_le : 0;
		//byte_order_is_be = ((uint8_t*)&magic_number)[i] == (uint8_t)((uint64_t)0xFFu & (magic_number >> (uint64_t)(((sizeof(magic_number)-(size_t)1u)-i)*(size_t)8u))) ? byte_order_is_be : 0;
	}

	if( (!byte_order_is_le && !byte_order_is_be) || (byte_order_is_be && byte_order_is_le) ) {
		fprintf(stderr, "%s error: Unable to detect host byte order (non-standard endianess).\n", "ply_write()");
		return STARFLOOD_FAILURE;
	}

	fprintf(ply_file,
		"ply\n"
		"format %s 1.0\n"
		"comment Created by Starflood version %d.%d.%d\n"
		"element vertex %zu\n"
		"property float x\n"
		"property float y\n"
		"property float z\n"
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
