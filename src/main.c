#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#pragma message("OpenMP is ENABLED")
#include <omp.h>
#else
#pragma message("OpenMP is NOT ENABLED")
#endif

#include "simulation.h"
#include "types.h"
#include "visualization.h"

int main(void) {
	printf("=== Starflood ===\n");

	fflush(stdout);

	printf("Configuration:\n");

	#ifdef _OPENMP
	printf("  OpenMP is ON (_OPENMP = %d)\n", _OPENMP);
	#else
	printf("  OpenMP is OFF\n");
	#endif

	printf("\n");

	unsigned int num_timesteps = 450u;

	//printf("");

	int rendering_enabled = true;

	simulation_t sim;

	sim.mem = NULL;

	if(simulation_init(&sim, 65536u) != EXIT_SUCCESS) {
		fprintf(stderr, "Error starting simulation!\n");

		return EXIT_FAILURE;
	}

	char filename[64];

	//simulation_load(&sim, filename);

	visualization_t vis;

	if(rendering_enabled) {
		if(visualization_init(&vis, 960u, 540u) != EXIT_SUCCESS) {
			fprintf(stderr, "Error starting visualization!\n");

			simulation_free(&sim);

			return EXIT_FAILURE;
		}
	}

	for(unsigned int step_num = 0u; step_num < num_timesteps; step_num++) {
		printf("Step #%3u\n", step_num);

		/*
		snprintf(filename, (size_t)64u, "./out/step_%04u.data", step_num);

		if(EXIT_SUCCESS != simulation_dump(&sim, filename)) {
			fprintf(stderr, "Error: simulation_dump() failed!\n");
		}
		*/

		if(rendering_enabled) {
			if(EXIT_SUCCESS != visualization_draw(&vis, &sim)) {
				fprintf(stderr, "Error: visualization_draw() failed!\n");
			}

			snprintf(filename, (size_t)64u, "./out/step_%04u.pfm", step_num);

			if(EXIT_SUCCESS != visualization_save(&vis, filename)) {
				fprintf(stderr, "Error: visualization_save() failed!\n");
			}
		}

		simulation_step(&sim);

		printf("\n");
	}

	if(rendering_enabled) {
		visualization_free(&vis);
	}

	simulation_free(&sim);

	return EXIT_SUCCESS;
}
