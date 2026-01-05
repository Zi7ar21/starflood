#include "grid.h"

#include <stdio.h>
#include <stdlib.h>

#include "common.h"
//#include "timing.h"

int grid_init(grid_t* restrict grid_ptr, unsigned int side_length, const real* restrict bounds_min, const real* restrict bounds_max) {
	//TIMING_INIT();

	grid_t grid = *grid_ptr;

	void* mem = NULL;

	if(  1u > side_length) {
		return STARFLOOD_FAILURE;
	}

	if(256u < side_length) {
		return STARFLOOD_FAILURE;
	}

	printf("%u\n", side_length);
	fflush(stdout);

	for(int i = 0; i < 3; i++) {
		grid.dim[i] = side_length;
	}

	for(int i = 0; i < 3; i++) {
		grid.bounds_min[i] = bounds_min[i];
		grid.bounds_max[i] = bounds_max[i];
	}

	size_t mem_length = (size_t)2u * (size_t)grid.dim[2] * (size_t)grid.dim[1] * (size_t)grid.dim[0];

	mem = malloc(sizeof(real) * mem_length);

	if(NULL == grid.mem) {
		return STARFLOOD_FAILURE;
	}

	{
		grid.mem = mem;
	}

	*grid_ptr = grid;

	return STARFLOOD_SUCCESS;
}

int grid_free(grid_t* restrict grid_ptr) {
	grid_t grid = *grid_ptr;

	free(grid.mem);

	for(int i = 0; i < 3; i++) {
		grid.dim[i] = 0u;
	}

	grid.mem = NULL;

	*grid_ptr = grid;

	return STARFLOOD_SUCCESS;
}

real* grid_samp(const grid_t* restrict grid_ptr) {
	grid_t grid = *grid_ptr;

	return (real*)grid.mem;
}
