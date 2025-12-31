#include "ply.h"

#include <stddef.h>
#include <stdio.h>

#include "common.h"

int ply_write(const char* restrict filename, float* restrict vertex_xyz, size_t num_vertices) {
	FILE* ply_file = fopen(filename, "wb");

	if(NULL == (void*)ply_file) {
		return STARFLOOD_FAILURE;
	}

	fprintf(ply_file,
		"ply\n"
		"format binary_little_endian 1.0\n"
		"comment Created by Starflood %d.%d.%d\n"
		"element vertex %zu\n"
		"property float x\n"
		"property float y\n"
		"property float z\n"
		"end_header\n",
		STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH,
		num_vertices
	);

	fflush(ply_file);

	fwrite(vertex_xyz, sizeof(float), (size_t)3u * num_vertices, ply_file);

	if( 0 != fclose(ply_file) ) {
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}
