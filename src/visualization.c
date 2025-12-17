// Needed for posix_memalign()
#define _POSIX_C_SOURCE 200112L

#include "visualization.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "config.h"
#include "rng.h"
#include "simulation.h"
#include "types.h"

int visualization_init(visualization_t* visualization, unsigned int w, unsigned int h) {
	visualization_t vis = *visualization;

	void* mem = NULL;

	i32* atomic_buffer = (i32*)NULL;
	f32* render_buffer = (f32*)NULL;

	if( sizeof(i32) != sizeof(f32) ) {
		return EXIT_FAILURE;
	}

	printf("Visualization Memory Addresses:\n");

	size_t atomic_buffer_length = (size_t)4u * (size_t)w * (size_t)h;
	size_t render_buffer_length = (size_t)4u * (size_t)w * (size_t)h;

	size_t atomic_buffer_offset = (size_t)0u;
	size_t render_buffer_offset = atomic_buffer_offset + atomic_buffer_length;

	size_t mem_size = (sizeof(i32)*atomic_buffer_length)+(sizeof(i32)*render_buffer_length);

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
	#endif

	if(NULL == mem) {
		#ifdef STARFLOOD_ALIGNMENT
		fprintf(stderr, "error in posix_memalign(&mem, %zu, %zu) while allocating memory for the ", (size_t)STARFLOOD_ALIGNMENT, mem_size);
		#else
		fprintf(stderr, "error in malloc(%zu) while allocating memory for the ", mem_size);
		#endif

		perror("visualization");

		return EXIT_FAILURE;
	}

	printf("  mem: %p\n", mem);

	memset(mem, 0, mem_size);

	atomic_buffer = (i32*)mem + (size_t)atomic_buffer_offset;
	render_buffer = (f32*)mem + (size_t)render_buffer_offset;

	printf("  atomic_buffer: %p (+%zu)\n", (void*)atomic_buffer, atomic_buffer_offset);
	printf("  render_buffer: %p (+%zu)\n", (void*)render_buffer, render_buffer_offset);
	printf("\n");

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = (i32)0;
	}

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		render_buffer[i] = (f32)0.000;
	}

	vis.w = w;
	vis.h = h;
	vis.mem = mem;
	vis.atomic_buffer = atomic_buffer;
	vis.render_buffer = render_buffer;

	*visualization = vis;

	return EXIT_SUCCESS;
}

int visualization_draw(visualization_t* visualization, simulation_t* simulation) {
	#ifdef _OPENMP
	double t0 = omp_get_wtime();
	double t1 = omp_get_wtime();
	#endif

	simulation_t sim = *simulation;

	visualization_t vis = *visualization;

	unsigned int w = vis.w;
	unsigned int h = vis.h;

	i32* atomic_buffer = vis.atomic_buffer;
	f32* render_buffer = vis.render_buffer;

	unsigned int N = sim.N;

	unsigned int step_number = sim.step_number;

	real* pos = sim.pos;

	#ifdef _OPENMP
	t0 = omp_get_wtime();
	#endif

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = (i32)0;
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("visualization_draw: %.09f ms clear atomic_buffer\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	/*
	double tran[16u] = {
		 1.000,  0.000,  0.000,  tr[0],
		 0.000,  1.000,  0.000,  tr[1],
		 0.000,  0.000,  1.000,  tr[2],
		 0.000,  0.000,  0.000,  1.000
	};

	double view[16u] = {
		 1.000,  0.000,  0.000,  0.000,
		 0.000,  1.000,  0.000,  0.000,
		 0.000,  0.000,  1.000,  0.000,
		 0.000,  0.000,  0.000,  1.000
	};

	real matrix[16];
	*/

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic, 1024)
	#endif
	for(unsigned int id = 0u; id < N; id++) {
		// https://cs418.cs.illinois.edu/website/text/math2.html
		double p[4] = {
			(double)pos[3u*id+0u],
			(double)pos[3u*id+1u],
			(double)pos[3u*id+2u],
			(double)1.0
		};

		/*
		{
			real mv[4] = {
				(matrix[ 0]*vector[0])+(matrix[ 1]*vector[1])+(matrix[ 2]*vector[2])+(matrix[ 3]*vector[3]),
				(matrix[ 4]*vector[0])+(matrix[ 5]*vector[1])+(matrix[ 6]*vector[2])+(matrix[ 7]*vector[3]),
				(matrix[ 8]*vector[0])+(matrix[ 9]*vector[1])+(matrix[10]*vector[2])+(matrix[11]*vector[3]),
				(matrix[12]*vector[0])+(matrix[13]*vector[1])+(matrix[14]*vector[2])+(matrix[15]*vector[3]),
			};

			p[0] = p_new[0];
			p[1] = p_new[1];
			p[2] = p_new[2];
			p[3] = p_new[3];
		}
		*/

		for(unsigned int j = 0u; j < 4u; j++) {
			p[j] *= exp2(ORTHO_SCALE);
		}

		// rotate2d
		{
			//const double theta = 0.0000 * TAU;
			//const double theta = 0.0625 * TAU;
			//const double theta = 0.1250 * TAU;
			const double theta = 0.2500 * TAU;

			double v[2] = {
				p[1], // Y
				p[2]  // Z
			};

			double vr[2] = {
				cos(theta)*v[0]-sin(theta)*v[1],
				sin(theta)*v[0]+cos(theta)*v[1]
			};

			p[1] = vr[0];
			p[2] = vr[1];
		}

		i32 rgba[4] = {(i32)1, (i32)1, (i32)1, (i32)1};

		for(unsigned int j = 0u; j < (unsigned int)SPATIAL_SAMPLES; j++) {
			u32 s[4] = {(u32)j, (u32)id, (u32)step_number, (u32)0xCF52CA01u};

			pcg4d(s);
			//pcg4d(s);

			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Gaussian pixel filter
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			double sample_offset[2] = {
				0.750 * sqrt( -2.0 * log(r[0]) ) * cos(TAU * r[1]),
				0.750 * sqrt( -2.0 * log(r[0]) ) * sin(TAU * r[1])
			};

			int x = (int)( ( (double)h * (0.500*p[0]) )+(0.5*(double)w)+sample_offset[0]-0.5 );
			int y = (int)( ( (double)h * (0.500*p[1]) )+(0.5*(double)h)+sample_offset[1]-0.5 );

			if((int)0 <= x && x < (int)w && (int)0 <= y && y < (int)h) {
				size_t pixel_index = (size_t)4u*((size_t)w*(size_t)y+(size_t)x);

				for(size_t k = (size_t)0u; k < (size_t)4u; k++) {
					#ifdef _OPENMP
					#pragma omp atomic update
					#endif
					atomic_buffer[pixel_index + k] += rgba[k];
				}
			}
		}
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("visualization_draw: %.09f ms rasterization\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	// Since atomic_buffer colors are quantized integers, it needs a floating-point scaling factor
	double pixel_value_scale = 1.0;

	#ifdef SPATIAL_SAMPLES
	pixel_value_scale *= 1.0 / (double)SPATIAL_SAMPLES;
	#endif

	#ifdef EXPOSURE
	pixel_value_scale *= exp2(EXPOSURE);
	#endif

	#ifdef _OPENMP
	t0 = omp_get_wtime();
	#endif

	// Read the atomic buffer and finish rendering the visualization
	for(unsigned int i = 0u; i < w * h; i++) {
		double rgba[4] = {
			(double)0.250,
			(double)0.250,
			(double)0.250,
			(double)1.000
		};

		for(unsigned int j = 0u; j < 4u; j++) {
			i32 val = (i32)0;

			#ifdef _OPENMP
			#pragma omp atomic read
			#endif
			val = atomic_buffer[4u*i+j];

			rgba[j] = pixel_value_scale * (double)val;
		}

		// Tonemapping
		for(unsigned int j = 0u; j < 4u; j++) {
			rgba[j] = tanh(rgba[j]);
		}

		// Set alpha channel to 1.000
		rgba[3u] = (double)1.000;

		for(unsigned int j = 0u; j < 4u; j++) {
			render_buffer[4u*i+j] = (f32)rgba[j];
		}
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("visualization_draw: %.09f ms post-processing\n", 1000.0*(t1-t0));
	#endif

	*visualization = vis;

	return EXIT_SUCCESS;
}

int visualization_save(visualization_t* visualization, const char* restrict filename) {
	#ifdef _OPENMP
	double t0 = omp_get_wtime();
	double t1 = omp_get_wtime();
	#endif

	visualization_t vis = *visualization;

	unsigned int image_w = vis.w;
	unsigned int image_h = vis.h;

	f32* render_buffer = vis.render_buffer;

	#ifdef _OPENMP
	t0 = omp_get_wtime();
	#endif

	FILE* file = fopen(filename, "wb");

	if(NULL == file) {
		fprintf(stderr, "Error: fopen(%s, \"wb\") failed: ", filename);

		perror("Error");

		return EXIT_FAILURE;
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("visualization_save: %.09f ms fopen()\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	// PFM graphic image file format
	// https://netpbm.sourceforge.net/doc/pfm.html
	fprintf(file, "PF\n%zu %zu\n-1.0\n", (size_t)image_w, (size_t)image_h);

	for(unsigned int y = 0u; y < image_h; y++) {
		for(unsigned int x = 0u; x < image_w; x++) {
			f32 rgb[3] = {
				render_buffer[4u*(image_w*y+x)+0u],
				render_buffer[4u*(image_w*y+x)+1u],
				render_buffer[4u*(image_w*y+x)+2u]
			};

			fwrite(rgb, sizeof(f32), (size_t)3u, file);
		}
	}

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("visualization_save: %.09f ms fwrite() loop\n", 1000.0*(t1-t0));

	t0 = omp_get_wtime();
	#endif

	fclose(file);

	#ifdef _OPENMP
	t1 = omp_get_wtime();

	printf("visualization_save: %.09f ms fclose()\n", 1000.0*(t1-t0));
	#endif

	return EXIT_SUCCESS;
}

int visualization_free(visualization_t* visualization) {
	visualization_t vis = *visualization;

	{
		void* mem = vis.mem;

		free(mem);
	}

	vis.mem = NULL;

	vis.atomic_buffer = (i32*)NULL;
	vis.render_buffer = (f32*)NULL;

	*visualization = vis;

	return EXIT_SUCCESS;
}
