// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "grid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "timing.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_FFT
#include <fftw3.h>
#endif

int grid_init(grid_t* restrict grid_ptr, unsigned int side_length, const real* restrict bounds_min, const real* restrict bounds_max) {
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

	if(NULL == mem) {
		return STARFLOOD_FAILURE;
	}

	#ifdef ENABLE_FFT
		#ifdef STARFLOOD_DOUBLE_PRECISION
		//fftw_init_threads();
		fftw_import_wisdom_from_filename(OUTPUT_DIR "/fftw_wisdom");
		void* mem_fft0 = fftw_malloc(sizeof(fftw_complex) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2]);
		#else
		//fftwf_init_threads();
		fftwf_import_wisdom_from_filename(OUTPUT_DIR "/fftw_wisdom");
		void* mem_fft0 = fftwf_malloc(sizeof(fftwf_complex) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2]);
		#endif

		if(NULL == mem_fft0) {
			free(mem);
			return STARFLOOD_FAILURE;
		}

		#ifdef STARFLOOD_DOUBLE_PRECISION
		void* mem_fft1 = fftw_malloc(sizeof(fftw_complex) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2]);
		void* mem_fft2 = fftw_malloc(sizeof(fftw_complex) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2]);
		#else
		void* mem_fft1 = fftwf_malloc(sizeof(fftwf_complex) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2]);
		void* mem_fft2 = fftwf_malloc(sizeof(fftwf_complex) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2]);
		#endif
	#endif

	{
		grid.mem     = mem;
		#ifdef ENABLE_FFT
		grid.mem_fft0 = mem_fft0;
		grid.mem_fft1 = mem_fft1;
		grid.mem_fft2 = mem_fft2;
		#endif
	}

	*grid_ptr = grid;

	return STARFLOOD_SUCCESS;
}

int grid_free(grid_t* restrict grid_ptr) {
	grid_t grid = *grid_ptr;

	free(grid.mem);

	#ifdef ENABLE_FFT
		#ifdef STARFLOOD_DOUBLE_PRECISION
		fftw_free(grid.mem_fft2);
		fftw_free(grid.mem_fft1);
		fftw_free(grid.mem_fft0);
		//fftw_cleanup_threads();
		fftw_export_wisdom_to_filename(OUTPUT_DIR "/fftw_wisdom");
		#else
		fftwf_free(grid.mem_fft2);
		fftwf_free(grid.mem_fft1);
		fftwf_free(grid.mem_fft0);
		//fftwf_cleanup_threads();
		fftwf_export_wisdom_to_filename(OUTPUT_DIR "/fftw_wisdom");
		#endif
	#endif

	for(int i = 0; i < 3; i++) {
		grid.dim[i] = 0u;
	}

	grid.mem = NULL;

	#ifdef ENABLE_FFT
	grid.mem_fft0 = NULL;
	grid.mem_fft1 = NULL;
	grid.mem_fft2 = NULL;
	#endif

	*grid_ptr = grid;

	return STARFLOOD_SUCCESS;
}

real* grid_samp(const grid_t* restrict grid_ptr) {
	grid_t grid = *grid_ptr;

	return (real*)grid.mem;
}

#ifdef ENABLE_FFT
int grid_solve_fft(const grid_t* restrict grid_ptr) {
	TIMING_INIT();

	grid_t grid = *grid_ptr;

	unsigned int grid_dim[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

	real bounds_min[3] = {grid.bounds_min[0], grid.bounds_min[1], grid.bounds_min[2]};
	real bounds_max[3] = {grid.bounds_max[0], grid.bounds_max[1], grid.bounds_max[2]};

	real bounds_wid[3] = {
		bounds_max[0] - bounds_min[0],
		bounds_max[1] - bounds_min[1],
		bounds_max[2] - bounds_min[2]
	};

	real* samp = grid_samp(&grid);

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_plan p;

	fftw_complex* grid_in = (fftw_complex*)grid.mem_fft0;
	fftw_complex* fourier = (fftw_complex*)grid.mem_fft1;
	fftw_complex* stencil = (fftw_complex*)grid.mem_fft2;
	#else
	fftwf_plan p;

	fftwf_complex* grid_in = (fftwf_complex*)grid.mem_fft0;
	fftwf_complex* fourier = (fftwf_complex*)grid.mem_fft1;
	fftwf_complex* stencil = (fftwf_complex*)grid.mem_fft2;
	#endif

	TIMING_START();

	for(unsigned int i = 0u; i < grid_dim[2] * grid_dim[1] * grid_dim[0]; i++) {
		grid_in[i][0] = (real)0.0;
		grid_in[i][1] = (real)0.0;
	}

	for(unsigned int i = 0u; i < grid_dim[2] * grid_dim[1] * grid_dim[0]; i++) {
		fourier[i][0] = (real)0.0;
		fourier[i][1] = (real)0.0;
	}

	for(unsigned int i = 0u; i < grid_dim[2] * grid_dim[1] * grid_dim[0]; i++) {
		stencil[i][0] = (real)0.0;
		stencil[i][1] = (real)0.0;
	}

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "clear_mem_fft");

	real bounds_inv[3] = {
		(real)(1.0 / (double)bounds_wid[0]),
		(real)(1.0 / (double)bounds_wid[1]),
		(real)(1.0 / (double)bounds_wid[2])
	};

	TIMING_START();

	#ifdef _OPENMP
	//#pragma omp for collapse(3)
	#endif
	for(unsigned int z = 0u; z < grid.dim[2]; z++) {
		for(unsigned int y = 0u; y < grid.dim[1]; y++) {
			for(unsigned int x = 0u; x < grid.dim[0]; x++) {
				unsigned int index = (grid.dim[1] * grid.dim[0] * z) + (grid.dim[0] * y) + x;

				real pos[3] = {
					(bounds_wid[0] * ( (real)x / (real)grid_dim[0] ) + bounds_min[0] ),
					(bounds_wid[1] * ( (real)y / (real)grid_dim[1] ) + bounds_min[1] ),
					(bounds_wid[2] * ( (real)z / (real)grid_dim[2] ) + bounds_min[2] )
				};

				real r2 = (pos[0]*pos[0])+(pos[1]*pos[1])+(pos[2]*pos[2]);

				real inv_r1 = (real)0.0 < r2 ? (real)1.0 / real_sqrt(r2) : (real)0.0;

				real phi = (real)(-G) * inv_r1;

				real c[2] = {
					(real)phi,
					(real)0.0
				};

				grid_in[index][0] = c[0];
				grid_in[index][1] = c[1];
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "calculate_stencil");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	//#ifdef _OPENMP
	//fftw_plan_with_nthreads(omp_get_max_threads());
	//#endif
	p =  fftw_plan_dft_3d(grid.dim[0], grid.dim[1], grid.dim[2], grid_in, stencil, FFTW_FORWARD, FFTW_MEASURE);
	#else
	//#ifdef _OPENMP
	//fftwf_plan_with_nthreads(omp_get_max_threads());
	//#endif
	p = fftwf_plan_dft_3d(grid.dim[0], grid.dim[1], grid.dim[2], grid_in, stencil, FFTW_FORWARD, FFTW_MEASURE);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "stencil_fftw_plan_dft_3d");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_execute(p);
	#else
	fftwf_execute(p);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "stenctil_fftw_execute_plan");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_destroy_plan(p);
	#else
	fftwf_destroy_plan(p);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "stencil_fftw_destroy_plan");
	TIMING_START();

	#ifdef _OPENMP
	//#pragma omp for collapse(3)
	#endif
	for(unsigned int z = 0u; z < grid.dim[2]; z++) {
		for(unsigned int y = 0u; y < grid.dim[1]; y++) {
			for(unsigned int x = 0u; x < grid.dim[0]; x++) {
				unsigned int index = (grid.dim[1] * grid.dim[0] * z) + (grid.dim[0] * y) + x;

				real phi = samp[index];

				real c[2] = {
					(real)phi,
					(real)0.0
				};

				grid_in[index][0] = c[0];
				grid_in[index][1] = c[1];
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "grid_in");

	for(unsigned int i = 0u; i < grid_dim[2] * grid_dim[1] * grid_dim[0]; i++) {
		fourier[i][0] = (real)0.0;
		fourier[i][1] = (real)0.0;
	}

	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	//#ifdef _OPENMP
	//fftw_plan_with_nthreads(omp_get_max_threads());
	//#endif
	p =  fftw_plan_dft_3d(grid.dim[0], grid.dim[1], grid.dim[2], grid_in, fourier, FFTW_FORWARD, FFTW_MEASURE);
	#else
	//#ifdef _OPENMP
	//fftwf_plan_with_nthreads(omp_get_max_threads());
	//#endif
	p = fftwf_plan_dft_3d(grid.dim[0], grid.dim[1], grid.dim[2], grid_in, fourier, FFTW_FORWARD, FFTW_MEASURE);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "grid_in_fftw_plan_dft_3d");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_execute(p);
	#else
	fftwf_execute(p);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "grid_in_fftw_execute_plan");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_destroy_plan(p);
	#else
	fftwf_destroy_plan(p);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "grid_in_fftw_destroy_plan");
	TIMING_START();

	#ifdef _OPENMP
	//#pragma omp for collapse(3)
	#endif
	for(unsigned int z = 0u; z < grid.dim[2]; z++) {
		for(unsigned int y = 0u; y < grid.dim[1]; y++) {
			for(unsigned int x = 0u; x < grid.dim[0]; x++) {
				unsigned int index = (grid.dim[1] * grid.dim[0] * z) + (grid.dim[0] * y) + x;

				//real phi = samp[index];

				real c_fourier[2] = {
					fourier[index][0],
					fourier[index][1]
				};

				real c_stencil[2] = {
					stencil[index][0],
					stencil[index][1]
				};

				real c[2] = {
					(c_fourier[0]*c_stencil[0])-(c_fourier[1]*c_stencil[1]),
					(c_fourier[0]*c_stencil[1])+(c_fourier[1]*c_stencil[0])
				};

				fourier[index][0] = c[0];
				fourier[index][1] = c[1];
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "multiply_conjugate");

	for(unsigned int i = 0u; i < grid_dim[2] * grid_dim[1] * grid_dim[0]; i++) {
		grid_in[i][0] = (real)0.0;
		grid_in[i][1] = (real)0.0;
	}

	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	//#ifdef _OPENMP
	//fftw_plan_with_nthreads(omp_get_max_threads());
	//#endif
	p =  fftw_plan_dft_3d(grid.dim[0], grid.dim[1], grid.dim[2], fourier, grid_in, FFTW_BACKWARD, FFTW_MEASURE);
	#else
	//#ifdef _OPENMP
	//fftwf_plan_with_nthreads(omp_get_max_threads());
	//#endif
	p = fftwf_plan_dft_3d(grid.dim[0], grid.dim[1], grid.dim[2], fourier, grid_in, FFTW_BACKWARD, FFTW_MEASURE);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "fourier_fftw_plan_dft_3d");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_execute(p);
	#else
	fftwf_execute(p);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "fourier_fftw_execute_plan");
	TIMING_START();

	#ifdef STARFLOOD_DOUBLE_PRECISION
	fftw_destroy_plan(p);
	#else
	fftwf_destroy_plan(p);
	#endif

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "fourier_fftw_destroy_plan");
	TIMING_START();

	#ifdef _OPENMP
	//#pragma omp for collapse(3)
	#endif
	for(unsigned int z = 0u; z < grid.dim[2]; z++) {
		for(unsigned int y = 0u; y < grid.dim[1]; y++) {
			for(unsigned int x = 0u; x < grid.dim[0]; x++) {
				unsigned int coord[3] = {x, y, z};

				coord[0] = ( coord[0] + (grid.dim[0] / 2u) ) % grid.dim[0];
				coord[1] = ( coord[1] + (grid.dim[1] / 2u) ) % grid.dim[1];
				coord[2] = ( coord[2] + (grid.dim[2] / 2u) ) % grid.dim[2];

				unsigned int index = (grid.dim[1] * grid.dim[0] * coord[2]) + (grid.dim[0] * coord[1]) + coord[0];

				real c[2] = {
					(real)grid_in[index][0],
					(real)grid_in[index][1]
				};

				real c2 = (c[0] * c[0]) + (c[1] * c[1]);

				index = (grid.dim[1] * grid.dim[0] * z) + (grid.dim[0] * y) + x;

				//samp[index] = real_sqrt(c2);
				samp[index] = c[0];
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("grid_solve_fft()", "finish");
	TIMING_START();

	return STARFLOOD_SUCCESS;
}
#endif
