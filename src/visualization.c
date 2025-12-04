#include "visualization.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "rng.h"
#include "simulation.h"
#include "types.h"

int visualization_init(visualization_t* visualization, unsigned int w, unsigned int h) {
	visualization_t vis = *visualization;

	printf("visualization memory addresses:\n");

	size_t atomic_buffer_size = sizeof(i32)*(size_t)4u*(size_t)w*(size_t)h;

	i32* atomic_buffer = (i32*)aligned_alloc(STARFLOOD_ALIGNMENT, atomic_buffer_size);

	if(NULL == (void*)atomic_buffer) {
		fprintf(stderr, "error in aligned_alloc(%zu, %zu) while allocating atomic_buffer", STARFLOOD_ALIGNMENT, atomic_buffer_size);

		perror("");

		return EXIT_FAILURE;
	}

	printf("  atomic_buffer: %p\n", (void*)atomic_buffer);

	size_t render_buffer_size = sizeof(f32)*(size_t)4u*(size_t)w*(size_t)h;

	f32* render_buffer = (f32*)aligned_alloc(STARFLOOD_ALIGNMENT, render_buffer_size);

	if(NULL == (void*)render_buffer) {
		fprintf(stderr, "error in aligned_alloc(%zu, %zu) while allocating render_buffer", STARFLOOD_ALIGNMENT, atomic_buffer_size);

		perror("");

		free(atomic_buffer);

		return EXIT_FAILURE;
	}

	printf("  render_buffer: %p\n", (void*)atomic_buffer);
	printf("\n");

	memset((void*)atomic_buffer, 0, atomic_buffer_size);
	memset((void*)render_buffer, 0, render_buffer_size);

	vis.w = w;
	vis.h = h;
	vis.atomic_buffer = atomic_buffer;
	vis.render_buffer = render_buffer;

	*visualization = vis;

	return EXIT_SUCCESS;
}

int visualization_draw(visualization_t* visualization, simulation_t* simulation) {
	simulation_t sim = *simulation;

	visualization_t vis = *visualization;

	unsigned int w = vis.w;
	unsigned int h = vis.h;

	i32* atomic_buffer = vis.atomic_buffer;
	f32* render_buffer = vis.render_buffer;

	unsigned int N = sim.N;

	unsigned int step_number = sim.step_number;

	real* pos = sim.pos;

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = (i32)0;
	}

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic, 1024)
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		// https://cs418.cs.illinois.edu/website/text/math2.html
		real v[4] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u],
			(real)1.0
		};

		// Pitch
		{
			real theta = (real)TAU * 0.125;

			real v_rot[3] = {
				v[0],
				v[1]*cosf(theta)-v[2]*sinf(theta),
				v[1]*sinf(theta)+v[2]*cosf(theta),
			};

			v[0] = v_rot[0];
			v[1] = v_rot[1];
			v[2] = v_rot[2];
		}

		i32 rgba[4] = {(f32)1.000, (f32)1.000, (f32)1.000, (f32)1.000};

		for(unsigned int j = 0u; j < 64u; j++) {
			u32 s[4] = {(u32)j, (u32)i, (u32)step_number, (u32)420691337u};

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
				0.66666666666666666666666666666667 * sqrt( -2.0 * log(r[0]) ) * cos(TAU * r[1]),
				0.66666666666666666666666666666667 * sqrt( -2.0 * log(r[0]) ) * sin(TAU * r[1])
			};

			int x = ( h * (0.250 * (double)v[0]) )+(0.5 * (double)w) - 0.5 + sample_offset[0];
			int y = ( h * (0.250 * (double)v[1]) )+(0.5 * (double)h) - 0.5 + sample_offset[1];

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
	#pragma omp parallel for schedule(dynamic, 1024)
	#endif
	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		i32 val = (i32)0;

		#ifdef _OPENMP
		#pragma omp atomic read
		#endif
		val = atomic_buffer[i];

		render_buffer[i] = (f32)tanh(0.125 * 0.125 * (double)val);
	}

	*visualization = vis;

	return EXIT_SUCCESS;
}

int visualization_save(visualization_t* visualization, const char* restrict filename) {
	visualization_t vis = *visualization;

	unsigned int image_w = vis.w;
	unsigned int image_h = vis.h;

	f32* render_buffer = vis.render_buffer;

	FILE* file = fopen(filename, "wb");

	if(NULL == file) {
		fprintf(stderr, "Error: fopen(%s, \"wb\") failed: ", filename);

		perror("Error");

		return EXIT_FAILURE;
	}

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

	fclose(file);

	return EXIT_SUCCESS;
}

int visualization_free(visualization_t* visualization) {
	visualization_t vis = *visualization;

	free(vis.render_buffer);
	free(vis.atomic_buffer);

	vis.atomic_buffer = NULL;
	vis.render_buffer = NULL;

	*visualization = vis;

	return EXIT_SUCCESS;
}
