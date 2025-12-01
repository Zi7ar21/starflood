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

#include "simulation.h"
#include "types.h"
#include "visualization.h"

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

	printf("\n");

	unsigned int num_timesteps = 10u;

	//printf("");

	int rendering_enabled = true;

	unsigned int N = 20000u;

	simulation_t sim;

	if(simulation_init(&sim, N) != EXIT_SUCCESS) {
		fprintf(stderr, "Error starting simulation!\n");

		return EXIT_FAILURE;
	}

	visualization_t vis;

	if(rendering_enabled) {
		if(visualization_init(&vis, 480u, 480u) != EXIT_SUCCESS) {
			fprintf(stderr, "Error starting visualization!\n");

			simulation_free(&sim);

			return EXIT_FAILURE;
		}
	}

	for(unsigned int n = 0u; n < num_timesteps; n++) {
		if(rendering_enabled) {
			if(visualization_draw(&vis, &sim) != EXIT_SUCCESS) {
			}
		}

		simulation_step(&sim);
	}

	if(rendering_enabled) {
		visualization_free(&vis);
	}

	simulation_free(&sim);

	return EXIT_SUCCESS;
}
