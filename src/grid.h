#pragma once

#include "config.h"
#include "types.h"

enum grid_conf {
	GRID_POTENTIAL,
	GRID_DOF
};

typedef struct {
	unsigned int dim[3];

	real bounds_min[3];
	real bounds_max[3];

	void* mem;

	#ifdef ENABLE_FFT
	void* mem_fft0;
	void* mem_fft1;
	void* mem_fft2;
	#endif
} grid_t;

// Initializes a grid
int grid_init(grid_t* restrict grid_ptr, unsigned int side_length, const real* restrict bounds_min, const real* restrict bounds_max);

// Frees a grid
int grid_free(grid_t* restrict grid_ptr);

// Get a samplable version of a grid
real* grid_samp(const grid_t* restrict grid_ptr);

#ifdef ENABLE_FFT
int grid_solve_fft(const grid_t* restrict grid_ptr);
#endif
