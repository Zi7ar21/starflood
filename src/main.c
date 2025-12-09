#include <math.h>
#include <stdbool.h>
#include <stddef.h>
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

#include "config.h"
#include "simulation.h"
#include "types.h"
#include "visualization.h"

int main(int argc, char *argv[]) {
	unsigned int N = (unsigned int)NUM_BODIES;
	unsigned int num_timesteps = (unsigned int)NUM_TIMESTEPS;

	int dumping_enabled = 0, visualization_enabled = 0;

	#ifdef ENABLE_DUMPING
	dumping_enabled = 1;
	#endif

	#ifdef ENABLE_VISUALIZATION
	visualization_enabled = 1;
	#endif

	unsigned int render_w = 960u, render_h = 960u;

	for(int i = 0; i < argc; i++) {
		if(NULL == argv[i]) {
			fprintf(stderr, "%s: error: argv[%d] is NULL!\n", argv[0], i);

			return EXIT_FAILURE;
		}

		if(0 == i) {
			continue;
		}

		if( 0 == strcmp("-h", argv[i]) || 0 == strcmp("--help", argv[i]) ) {
			printf(
				"Usage: %s [options]...\n"
				"\n"
				"Options:\n"
				"  -h, --help       display this help and exit\n"
				"  -V, --version    display version information and exit\n"
				"\n"
				"Starflood repository: <%s>\n",
				argv[0],
				STARFLOOD_REPOSITORY
			);

			if( (argc - 1) == i ) {
				return EXIT_SUCCESS;
			} else {
				return EXIT_FAILURE;
			}
		}

		if( 0 == strcmp("-V", argv[i]) || 0 == strcmp("--version", argv[i]) ) {
			printf(
				"Starflood %d.%d.%d\n"
				"\n"
				"Starflood repository: <%s>\n",
				STARFLOOD_VERSION_MAJOR,
				STARFLOOD_VERSION_MINOR,
				STARFLOOD_VERSION_PATCH,
				STARFLOOD_REPOSITORY
			);

			if( (argc - 1) == i ) {
				return EXIT_SUCCESS;
			} else {
				return EXIT_FAILURE;
			}
		}

		fprintf(stderr, "%s: unrecognized option \'%s\'\nTry \'%s --help\' for more information.\n", argv[0], argv[1], argv[0]);

		return EXIT_FAILURE;
	}

	printf("=== Starflood ===\n");

	printf("  Version %d.%d.%d\n", STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH);

	fflush(stdout);

	printf("\n");
	printf("Configuration:\n");
	printf("    OpenMP Support:\n");
	printf("        Enabled: ");
	#ifdef _OPENMP
	printf("true\n");
	printf("        _OPENMP: %d\n", _OPENMP);
	printf("        max_threads: %d\n", omp_get_max_threads());
	printf("        num_devices: %d\n", omp_get_num_devices());
	#else
	printf("false\n");
	#endif

	printf("Parameters:\n");
	printf("    Simulation:\n");
	printf("        N: %u\n", N);
	printf("        num_timesteps: %u\n", num_timesteps);
	printf("    Visualization:\n");
	printf("        Enabled: ");

	if(visualization_enabled) {
		printf("true\n");
	} else {
		printf("false\n");
	}

	printf("        Size: %ux%u\n", render_w, render_h);
	printf("\n");

	fflush(stdout);

	simulation_t sim;

	sim.mem = NULL;

	if(simulation_init(&sim, N) != EXIT_SUCCESS) {
		fprintf(stderr, "Error starting simulation!\n");

		return EXIT_FAILURE;
	}

	char filename[64];

	//simulation_load(&sim, filename);

	visualization_t vis;

	if(visualization_enabled) {
		if( EXIT_SUCCESS != visualization_init(&vis, render_w, render_h) ) {
			fprintf(stderr, "Error starting visualization!\n");

			simulation_free(&sim);

			return EXIT_FAILURE;
		}
	}

	for(unsigned int step_num = 0u; step_num <= num_timesteps; step_num++) {
		printf("Step #%3u\n", step_num);

		if(dumping_enabled) {
			snprintf(filename, (size_t)64u, SIMULATION_FILENAME, step_num);

			if( EXIT_SUCCESS != simulation_dump(&sim, filename) ) {
				fprintf(stderr, "Error: simulation_dump() failed!\n");
			}
		}

		if(visualization_enabled) {
			if( EXIT_SUCCESS != visualization_draw(&vis, &sim) ) {
				fprintf(stderr, "Error: visualization_draw() failed!\n");
			}

			snprintf(filename, (size_t)64u, VISUALIZATION_FILENAME, step_num);

			if( EXIT_SUCCESS != visualization_save(&vis, filename) ) {
				fprintf(stderr, "Error: visualization_save() failed!\n");
			}
		}

		if(num_timesteps != step_num) {
			simulation_step(&sim);
		}

		printf("\n");
	}

	printf("Cleaning up...\n");

	fflush(stdout);

	if(visualization_enabled) {
		visualization_free(&vis);
	}

	simulation_free(&sim);

	printf("Finished!\n");

	return EXIT_SUCCESS;
}
