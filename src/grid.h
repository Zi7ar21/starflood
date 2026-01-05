#pragma once

#include "types.h"

typedef struct {
	unsigned int dim[3];

	real bounds_min[3];
	real bounds_max[3];

	void* mem;

	#ifdef ENABLE_FFT
	void* mem_fft;
	#endif
} grid_t;

int grid_init(grid_t* restrict grid_ptr, unsigned int side_length, const real* restrict bounds_min, const real* restrict bounds_max);

int grid_free(grid_t* restrict grid_ptr);
