// needed for clock_getres() and clock_gettime()
#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
//#pragma message("OpenMP is ENABLED")
#include <omp.h>
#else
//#pragma message("OpenMP is NOT ENABLED")
#endif

#include "common.h"
#include "config.h"
#include "simulation.h"
#include "types.h"
#include "timing.h"
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
	printf("        wtick: %.03f ns\n", 1.0e9*omp_get_wtick());
	#else
	printf("false\n");
	#endif
	printf("    Timing:\n");
	printf("        Enabled: ");
	#ifdef ENABLE_TIMING
		printf("true\n");
		printf("        POSIX clock enabled: ");
		#ifdef TIMING_USE_OMP_GET_WTIME
		printf("false\n");
		#else
		printf("true\n");
		printf("            clockid: ");
		switch(STARFLOOD_POSIX_CLOCKID) {
			case CLOCK_REALTIME:
				printf("CLOCK_REALTIME");
				break;
			case CLOCK_PROCESS_CPUTIME_ID:
				printf("CLOCK_PROCESS_CPUTIME_ID");
				break;
			case CLOCK_MONOTONIC_RAW:
				printf("CLOCK_MONOTONIC_RAW");
				break;
			default:
				printf("Other/Unknown (%d)", STARFLOOD_POSIX_CLOCKID);
				break;
		}
		printf("\n");
		printf("            res: ");
		{
			struct timespec res;

			if( 0 != clock_getres(STARFLOOD_POSIX_CLOCKID, &res) ) {
				printf("error");
			} else {
				printf("%jd nanoseconds", (intmax_t)1000000000l*(intmax_t)res.tv_sec+(intmax_t)res.tv_sec);
			}
		}
		printf("\n");
		printf("            instantaneously measured difference: ");
		{
			struct timespec tp0, tp1;

			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);

			printf("%jd nanoseconds", (intmax_t)1000000000l*((intmax_t)(tp1.tv_sec)-(intmax_t)(tp0.tv_sec))+((intmax_t)(tp1.tv_nsec)-(intmax_t)(tp0.tv_nsec)));
		}
		printf("\n");
		#endif
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
		if( STARFLOOD_SUCCESS != visualization_init(&vis, visualization_dimensions[0], visualization_dimensions[1]) ) {
			fprintf(stderr, "Error: %s failed!\n", "visualization_init()");

			return EXIT_FAILURE;
		}
	}

	if(enable_simulation || enable_simulation_io) {
		if( STARFLOOD_SUCCESS != simulation_init(&sim, num_bodies) ) {
			fprintf(stderr, "Error: %s failed!\n", "simulation_init()");

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
			#ifdef OUTPUT_INTERVAL
			snprintf(sim_filename, filename_size, SIMULATION_FILENAME, step_num / (unsigned int)OUTPUT_INTERVAL);
			#else
			snprintf(sim_filename, filename_size, SIMULATION_FILENAME, step_num);
			#endif
		#endif

		#ifdef VISUALIZATION_FILENAME
			#ifdef OUTPUT_INTERVAL
			snprintf(vis_filename, filename_size, VISUALIZATION_FILENAME, step_num / (unsigned int)OUTPUT_INTERVAL);
			#else
			snprintf(vis_filename, filename_size, VISUALIZATION_FILENAME, step_num);
			#endif
		#endif

		#ifdef OUTPUT_INTERVAL
		if( (0u == (step_num % OUTPUT_INTERVAL)) && enable_simulation_io ) {
		#else
		if(enable_simulation_io) {
		#endif
			if(enable_simulation) {
				if( STARFLOOD_SUCCESS != simulation_save(&sim, sim_filename) ) {
					fprintf(stderr, "Error: %s failed!\n", "simulation_save()");
				}
			} else {
				if( STARFLOOD_SUCCESS != simulation_read(&sim, sim_filename) ) {
					fprintf(stderr, "Error: %s failed!\n", "simulation_read()");
				}
			}
		}

		#ifdef OUTPUT_INTERVAL
		if( (0u == (step_num % OUTPUT_INTERVAL)) && enable_visualization) {
		#else
		if(enable_visualization) {
		#endif
			if( STARFLOOD_SUCCESS != visualization_draw(&vis, &sim) ) {
				fprintf(stderr, "Error: %s failed!\n", "visualization_draw()");
			}

			if(enable_visualization_io) {
				if( STARFLOOD_SUCCESS != visualization_save(&vis, vis_filename) ) {
					fprintf(stderr, "Error: %s failed!\n", "visualization_save()");
				}
			}
		}

		if(enable_simulation) {
			if( STARFLOOD_SUCCESS != simulation_step(&sim)) {
				fprintf(stderr, "Error: %s failed!\n", "simulation_step()");
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

	printf("Finished.\n");

	return EXIT_SUCCESS;
}
