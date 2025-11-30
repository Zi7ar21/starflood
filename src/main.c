#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _OPENMP
#pragma message("OpenMP is ENABLED")
#include <omp.h>
#else
#pragma message("OpenMP is NOT ENABLED")
#endif

#include "types.h"
#include "simulation.h"

const double TAU = 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341360;

int main(void) {
	printf("=== Starflood ===\n");

	fflush(stdout);

	printf("Configuration:\n");

	#ifdef _OPENMP
	printf("  OpenMP is ON (_OPENMP = %d)\n", _OPENMP);
	#else
	printf("  OpenMP is OFF\n");
	#endif

	unsigned int num_timesteps = 300u;

	//printf("");

	int rendering_enabled = true;

	unsigned int N = 20000u;

	simulation_t sim;

	if(simulation_init(&sim, N) != EXIT_SUCCESS) {
		fprintf(stderr, "Error creating simulation!");

		return EXIT_FAILURE;
	}

	size_t image_w = (size_t)480u;
	size_t image_h = (size_t)480u;

	float* image = (float*)aligned_alloc((size_t)4096u, sizeof(float)*(size_t)3*image_w*image_h);

	if(NULL == image) {
		perror("error allocating image");

		simulation_free(&sim);

		return EXIT_FAILURE;
	}

	printf("image: %p\n", (void*)image);

	for(size_t i = (size_t)0u; i < (size_t)3u*image_w*image_h; i++) {
		image[i] = (float)0.0;
	}

	for(unsigned int n = 0u; n < num_timesteps; n++) {
		if(rendering_enabled) {
			for(size_t y = (size_t)0u; y < image_h; y++) {
				for(size_t x = (size_t)0u; x < image_w; x++) {
					image[3*(image_w*y+x)+0] = (float)0;
					image[3*(image_w*y+x)+1] = (float)0;
					image[3*(image_w*y+x)+2] = (float)0;
				}
			}

			real* pos = sim.pos;

			for(unsigned int i = 0u; i < N; i++) {
				//double t = 0.125*0.125*sin(0.25*TAU*(1.0/15.0)*(double)n);

				//real u = pos[3u*i+0u]*cos(TAU*t)-pos[3u*i+2u]*sin(TAU*t);

				//int x = ( (float)image_h*((float)0.25*u           ) )+((float)0.5*(float)image_w)-(float)0.5;
				int x = ( (float)image_h*((float)0.25*pos[3u*i+0u]) )+((float)0.5*(float)image_w)-(float)0.5;
				int y = ( (float)image_h*((float)0.25*pos[3u*i+1u]) )+((float)0.5*(float)image_h)-(float)0.5;

				if((int)0 <= x && x < (int)image_w && (int)0 <= y && y < (int)image_h) {
					image[3*(image_w*y+x)+0] += (float)0.333;
					image[3*(image_w*y+x)+1] += (float)0.333;
					image[3*(image_w*y+x)+2] += (float)0.333;
				}
			}

			char filename[FILENAME_MAX];

			snprintf(filename, FILENAME_MAX, "./out/%04d.pfm", n);

			FILE* file = fopen(filename, "w");

			if(NULL == file) {
				fprintf(stderr, "fopen(%s, \"w\") failed", filename);

				perror("");

				break;
			}

			fprintf(file, "PF\n%zu %zu\n-1.0\n", image_w, image_h);

			for(size_t y = (size_t)0u; y < image_h; y++) {
				for(size_t x = (size_t)0u; x < image_w; x++) {
					float pixel[3] = {
						image[3*(image_w*y+x)+0],
						image[3*(image_w*y+x)+1],
						image[3*(image_w*y+x)+2]
					};

					fwrite(pixel, sizeof(float), (size_t)3, file);
				}
			}

			fclose(file);
		}

		simulation_step(&sim);
	}

	free(image);

	simulation_free(&sim);

	return EXIT_SUCCESS;
}
