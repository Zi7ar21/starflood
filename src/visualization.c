#include "visualization.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "types.h"
#include "simulation.h"

int visualization_init(visualization_t* visualization, unsigned int w, unsigned int h) {
	visualization_t vis = *visualization;

	printf("Visualization Memory Addresses:\n");

	size_t atomic_buffer_size = sizeof(int)*(size_t)3u*(size_t)w*(size_t)h;

	int* atomic_buffer = (int*)aligned_alloc(STARFLOOD_ALIGNMENT, atomic_buffer_size);

	if(NULL == (void*)atomic_buffer) {
		fprintf(stderr, "error in aligned_alloc(%zu, %zu) while allocating atomic_buffer", STARFLOOD_ALIGNMENT, atomic_buffer_size);

		perror("");

		return EXIT_FAILURE;
	}

	printf("  atomic_buffer: %p\n", (void*)atomic_buffer);

	size_t render_buffer_size = sizeof(float)*(size_t)3u*(size_t)w*(size_t)h;

	float* render_buffer = (float*)aligned_alloc(STARFLOOD_ALIGNMENT, render_buffer_size);

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

	vis.id = 0;
	vis.w  = w;
	vis.h  = h;
	vis.atomic_buffer = atomic_buffer;
	vis.render_buffer = render_buffer;

	*visualization = vis;

	return EXIT_SUCCESS;
}

int visualization_draw(visualization_t* visualization, simulation_t* simulation) {
	simulation_t sim = *simulation;

	visualization_t vis = *visualization;

	int id = vis.id;

	unsigned int w = vis.w;
	unsigned int h = vis.h;

	int* atomic_buffer = vis.atomic_buffer;

	float* render_buffer = vis.render_buffer;

	for(unsigned int i = 0u; i < 3u*w*h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = 0;
	}

	unsigned int N = sim.N;

	real* pos = sim.pos;

	for(unsigned int i = 0u; i < N; i++) {
		//double t = 0.125*0.125*sin(0.25*TAU*(1.0/15.0)*(double)n);

		//real u = pos[3u*i+0u]*cos(TAU*t)-pos[3u*i+2u]*sin(TAU*t);

		//int x = ( (float)image_h*((float)0.25*u           ) )+((float)0.5*(float)image_w)-(float)0.5;
		int x = ( (float)h*((float)0.25*(float)pos[3u*i+0u]) )+((float)0.5*(float)w)-(float)0.5;
		int y = ( (float)h*((float)0.25*(float)pos[3u*i+1u]) )+((float)0.5*(float)h)-(float)0.5;

		if((int)0 <= x && x < (int)w && (int)0 <= y && y < (int)h) {
			for(unsigned int j = 0u; j < 3u; j++) {
				size_t idx = (size_t)(3u*(w*(unsigned int)y+(unsigned int)x)+j);

				#ifdef _OPENMP
				#pragma omp atomic update
				#endif
				atomic_buffer[idx] += 1;
			}
		}
	}

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,1024)
	#endif
	for(unsigned int i = 0u; i < 3u*w*h; i++) {
		int val = 0;

		#ifdef _OPENMP
		#pragma omp atomic read
		#endif
		val = atomic_buffer[i];

		render_buffer[i] = (float)0.125*(float)val;
	}

	char filename[FILENAME_MAX];

	snprintf(filename, FILENAME_MAX, "./out/%04d.pfm", id);

	FILE* file = fopen(filename, "w");

	if(NULL == file) {
		fprintf(stderr, "fopen(%s, \"w\") failed", filename);

		perror("");

		return EXIT_FAILURE;
	}

	fprintf(file, "PF\n%zu %zu\n-1.0\n", (size_t)w, (size_t)h);

	for(unsigned int y = 0u; y < h; y++) {
		for(unsigned int x = 0u; x < w; x++) {
			float pixel[3] = {
				render_buffer[3u*(w*y+x)+0u],
				render_buffer[3u*(w*y+x)+1u],
				render_buffer[3u*(w*y+x)+2u]
			};

			fwrite(pixel, sizeof(float), (size_t)3u, file);
		}
	}

	fclose(file);

	id++;

	vis.id = id;

	*visualization = vis;

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
