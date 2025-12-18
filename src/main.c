#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
//#pragma message("OpenMP is ENABLED")
#include <omp.h>
#else
//#pragma message("OpenMP is NOT ENABLED")
#endif

#include "config.h"
#include "simulation.h"
#include "types.h"
#include "visualization.h"

int main(int argc, char** argv) {
	int enable_simulation = 1, enable_simulation_io = 0, enable_visualization = 0, enable_visualization_io = 0;

	#ifndef ENABLE_SIMULATION
	enable_simulation = 0;
	#endif

	#ifdef SIMULATION_FILENAME
	enable_simulation_io = 1;
	#endif

	#ifdef ENABLE_VISUALIZATION
	enable_visualization= 1;
	#endif

	#ifdef VISUALIZATION_FILENAME
	enable_visualization_io = 1;
	#endif

	unsigned int num_bodies = (unsigned int)NUM_BODIES, num_timesteps = (unsigned int)NUM_TIMESTEPS;

	unsigned int visualization_dimensions[2] = {1920u, 1080u};

	visualization_t vis;

	simulation_t sim;

	const size_t filename_size = 256;

	char sim_filename[256] = "", vis_filename[256] = "";

	#ifdef SIMULATION_FILENAME
	if( (filename_size - 1) > strlen(SIMULATION_FILENAME) ) {
		strcpy(sim_filename, SIMULATION_FILENAME);
	}
	#endif

	#ifdef VISUALIZATION_FILENAME
	if( (filename_size - 1) > strlen(VISUALIZATION_FILENAME) ) {
		strcpy(vis_filename, VISUALIZATION_FILENAME);
	}
	#endif

	for(int i = 0; i < argc; i++) {
		if(0 >= i) {
			continue;
		}

		if(NULL == (void*)argv[i]) {
			fprintf(stderr, "%s: error: argv[%d] is NULL!\n", argv[0], i);

			if(1 >= i) {
				break;
			}

			return EXIT_FAILURE;
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
			}

			continue;
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
			}

			continue;
		}

		fprintf(stderr, "%s: error: unrecognized option \'%s\'\nTry \'%s --help\' for more information.\n", argv[0], argv[1], argv[0]);

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
	printf("        wtick: %.09f ns\n", 1.0e9*omp_get_wtick());
	#else
	printf("false\n");
	#endif

	printf("Parameters:\n");
	printf("    Simulation:\n");
	printf("        enabled: %d\n", enable_simulation);
	printf("        I/O enabled: %d\n", enable_simulation_io);

	if(enable_simulation_io) {
		printf("        sim_filename: %s\n", sim_filename);
	}

	printf("        N: %u\n", num_bodies);

	#ifdef N_DIV
	printf("        N_DIV: %u\n", (unsigned int)N_DIV);
	#else
	printf("        N_DIV: %u\n", 0u);
	#endif

	printf("        num_timesteps: %u\n", num_timesteps);
	printf("    Visualization:\n");
	printf("        enabled: %d\n", enable_visualization);
	printf("        I/O enabled: %d\n", enable_visualization_io);

	if(enable_visualization_io) {
		printf("        vis_filename: %s\n", vis_filename);
	}

	printf("        dimensions: %ux%u\n", visualization_dimensions[0], visualization_dimensions[1]);
	printf("\n");

	fflush(stdout);

	if(enable_visualization) {
		if( EXIT_SUCCESS != visualization_init(&vis, visualization_dimensions[0], visualization_dimensions[1]) ) {
			fprintf(stderr, "Error starting visualization!\n");

			return EXIT_FAILURE;
		}
	}

	if(enable_simulation || enable_simulation_io) {
		if( EXIT_SUCCESS != simulation_init(&sim, num_bodies) ) {
			fprintf(stderr, "Error starting simulation!\n");

			if(enable_visualization) {
				visualization_free(&vis);
			}

			return EXIT_FAILURE;
		}
	}

	for(unsigned int step_num = 0u; step_num <= num_timesteps; step_num++) {
	//for(unsigned int step_num = 400u; step_num <= 400u; step_num++) {
		printf("Step #%3u\n", step_num);

		sim.step_number = step_num;

		#ifdef SIMULATION_FILENAME
		snprintf(sim_filename, filename_size, SIMULATION_FILENAME, step_num / (unsigned int)OUTPUT_INTERVAL);
		#endif

		#ifdef VISUALIZATION_FILENAME
		snprintf(vis_filename, filename_size, VISUALIZATION_FILENAME, step_num / (unsigned int)OUTPUT_INTERVAL);
		#endif

		if( (0u == (step_num % OUTPUT_INTERVAL)) && enable_simulation_io ) {
			if(enable_simulation) {
				if( EXIT_SUCCESS != simulation_save(&sim, sim_filename) ) {
					fprintf(stderr, "Error: simulation_save() failed!\n");
				}
			} else {
				if( EXIT_SUCCESS != simulation_read(&sim, sim_filename) ) {
					fprintf(stderr, "Error: simulation_read() failed!\n");
				}
			}
		}

		if( (0u == (step_num % OUTPUT_INTERVAL)) && enable_visualization) {
			if( EXIT_SUCCESS != visualization_draw(&vis, &sim) ) {
				fprintf(stderr, "Error: visualization_draw() failed!\n");
			}

			if(enable_visualization_io) {
				if( EXIT_SUCCESS != visualization_save(&vis, vis_filename) ) {
					fprintf(stderr, "Error: visualization_save() failed!\n");
				}
			}
		}

		if(enable_simulation) {
			if( EXIT_SUCCESS != simulation_step(&sim)) {
				fprintf(stderr, "Error: simulation_step() failed!\n");
			}
		}

		printf("\n");
	}

	printf("Cleaning up...\n");

	fflush(stdout);

	if(enable_simulation || enable_simulation_io) {
		simulation_free(&sim);
	}

	if(enable_visualization) {
		visualization_free(&vis);
	}

	printf("Finished!\n");

	return EXIT_SUCCESS;
}
