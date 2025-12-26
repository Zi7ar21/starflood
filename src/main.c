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


#define SIM_FILE_EXTENSION ".raw"

#if (0 == VISUALIZATION_IMAGE_FORMAT)
#define VIS_FILE_EXTENSION ".pfm"
#endif

#if (1 == VISUALIZATION_IMAGE_FORMAT)
#define VIS_FILE_EXTENSION ".ppm"
#endif

#if (2 == VISUALIZATION_IMAGE_FORMAT)
#define VIS_FILE_EXTENSION ".png"
#endif

int main(int argc, char** argv) {
	int enable_sim = 1, enable_sim_io = 0, enable_vis = 0, enable_vis_io = 0;

	#ifndef ENABLE_SIM
	enable_sim = 0;
	#endif

	#ifdef SIM_FILENAME
	enable_sim_io = 1;
	#endif

	#ifdef ENABLE_VIS
	enable_vis= 1;
	#endif

	#ifdef VIS_FILENAME
	enable_vis_io = 1;
	#endif

	unsigned int num_bodies = (unsigned int)NUM_BODIES, num_timesteps = (unsigned int)NUM_TIMESTEPS;

	// ( 1:1) Aspect Ratio
	//unsigned int visualization_dimensions[2] = { 144u,  144u}; // Potato (useful for computationally intensive phase-space searches)
	//unsigned int visualization_dimensions[2] = { 360u,  360u}; //   Standard Definition
	//unsigned int visualization_dimensions[2] = { 720u,  720u}; //       High-Definition
	//unsigned int visualization_dimensions[2] = {1080u, 1080u}; //  Full High-Definition
	//unsigned int visualization_dimensions[2] = {2160u, 2160u}; // Ultra High-Definition

	// (16:9) Aspect Ratio
	//unsigned int visualization_dimensions[2] = { 640u,  360u}; //   Standard Definition
	//unsigned int visualization_dimensions[2] = {1280u,  720u}; //       High-Definition
	unsigned int visualization_dimensions[2] = {1920u, 1080u}; //  Full High-Definition
	//unsigned int visualization_dimensions[2] = {3840u, 2160u}; // Ultra High-Definition

	// Cinematic Aspect Ratio
	//unsigned int visualization_dimensions[2] = {2048u, 1080u}; // DCI 2K
	//unsigned int visualization_dimensions[2] = {4096u, 2160u}; // DCI 4K

	sim_t sim;
	vis_t vis;

	const size_t filename_max = (size_t)STARFLOOD_FILENAME_MAX;

	/*
	char* output_dir = (char*)malloc(sizeof(char) * filename_max);

	if(NULL == (void*)output_dir) {
		fprintf(stderr, "fatal error: output_dir is (char*)NULL after malloc(%zu", sizeof(char) * filename_max);
		perror(")");
		return EXIT_FAILURE;
	}
	*/

	char* sim_filename_pattern = (char*)malloc(sizeof(char) * filename_max);

	if(NULL == (void*)sim_filename_pattern) {
		fprintf(stderr, "fatal error: sim_filename_pattern is (char*)NULL after malloc(%zu", sizeof(char) * filename_max);
		perror(")");
		return EXIT_FAILURE;
	}

	char* vis_filename_pattern = (char*)malloc(sizeof(char) * filename_max);

	if(NULL == (void*)vis_filename_pattern) {
		fprintf(stderr, "fatal error: vis_filename_pattern is (char*)NULL after malloc(%zu", sizeof(char) * filename_max);
		perror(")");
		free(sim_filename_pattern);
		return EXIT_FAILURE;
	}

	char* sim_filename = (char*)malloc(sizeof(char) * filename_max);

	if(NULL == (void*)sim_filename) {
		fprintf(stderr, "fatal error: sim_filename is (char*)NULL after malloc(%zu", sizeof(char) * filename_max);
		perror(")");
		free(sim_filename_pattern);
		free(vis_filename_pattern);
		return EXIT_FAILURE;
	}

	char* vis_filename = (char*)malloc(sizeof(char) * filename_max);

	if(NULL == (void*)vis_filename) {
		fprintf(stderr, "fatal error: vis_filename is (char*)NULL after malloc(%zu", sizeof(char) * filename_max);
		perror(")");
		free(sim_filename_pattern);
		free(vis_filename_pattern);
		free(sim_filename);
		return EXIT_FAILURE;
	}

	#ifdef SIM_FILENAME
	strcpy(sim_filename_pattern, OUTPUT_DIR "/" SIM_FILENAME SIM_FILE_EXTENSION);
	#else
	strcpy(sim_filename_pattern, OUTPUT_DIR "/sim/step_%04u" SIM_FILE_EXTENSION);
	#endif

	#ifdef VIS_FILENAME
	strcpy(vis_filename_pattern, OUTPUT_DIR "/" VIS_FILENAME VIS_FILE_EXTENSION);
	#else
	strcpy(vis_filename_pattern, OUTPUT_DIR "/vis/step_%04u" VIS_FILE_EXTENSION);
	#endif

	if( 0 >= snprintf(sim_filename, filename_max, sim_filename_pattern, 0u) ) {
		fprintf(stderr, "fatal error: snprintf(sim_filename, filename_max, \"%s\", %u) ", sim_filename_pattern, 0u);
		perror("failed");
		free(sim_filename);
		free(vis_filename);
		free(sim_filename_pattern);
		free(vis_filename_pattern);
		return EXIT_FAILURE;
	}

	if( 0 >= snprintf(vis_filename, filename_max, vis_filename_pattern, 0u) ) {
		fprintf(stderr, "fatal error: snprintf(vis_filename, filename_max, \"%s\", %u) ", vis_filename_pattern, 0u);
		perror("failed");
		free(sim_filename);
		free(vis_filename);
		free(sim_filename_pattern);
		free(vis_filename_pattern);
		return EXIT_FAILURE;
	}

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

	printf("\n  Version %d.%d.%d\n", STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH);

	fflush(stdout);

	printf("\nConfiguration:");
	printf("\n    OpenMP Support:");
	printf("\n        Enabled: ");
	#ifdef _OPENMP
	printf("true");
	printf("\n        _OPENMP: %d", _OPENMP);
	printf("\n        max_threads: %d", omp_get_max_threads());
	printf("\n        num_devices: %d", omp_get_num_devices());
	printf("\n        wtick: %.03f nanosecond(s)", 1.0e9*omp_get_wtick());
	#else
	printf("false");
	#endif
	printf("\n    Timing:");
	printf("\n        Enabled: ");
	#ifdef ENABLE_TIMING
		printf("true");
		printf("\n        POSIX clock enabled: ");
		#ifdef TIMING_USE_OMP_GET_WTIME
		printf("false");
		#else
		printf("true");
		printf("\n            clockid: ");

		switch(STARFLOOD_POSIX_CLOCKID) {
			case CLOCK_REALTIME:
				printf("CLOCK_REALTIME");
				break;
			case CLOCK_MONOTONIC:
				printf("CLOCK_MONOTONIC");
				break;
			case CLOCK_MONOTONIC_RAW:
				printf("CLOCK_MONOTONIC_RAW");
				break;
			default:
				printf("Other/Unknown (%d)", STARFLOOD_POSIX_CLOCKID);
				break;
		}

		{
			struct timespec tp0, tp1;

			printf("\n            res: ");

			if( 0 != clock_getres(STARFLOOD_POSIX_CLOCKID, &tp0) ) {
				printf("N/A\n");
			} else {
				printf("%jd nanosecond(s)", (intmax_t)1000000000l*(intmax_t)tp0.tv_sec+(intmax_t)tp0.tv_sec);
			}

			printf("\n            instantaneously measured difference: ");

			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp0);
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &tp1);

			printf("%jd nanosecond(s)", (intmax_t)1000000000l*((intmax_t)(tp1.tv_sec)-(intmax_t)(tp0.tv_sec))+((intmax_t)(tp1.tv_nsec)-(intmax_t)(tp0.tv_nsec)));

			printf("\n            current time: %jd.%09jd seconds", (intmax_t)tp0.tv_sec, (intmax_t)tp0.tv_nsec);
		}
		#endif
	#else
	printf("false");
	#endif
	printf("\nParameters:");
	printf("\n    num_timesteps: %u", num_timesteps);
	printf("\n    Simulation:");
	printf("\n        enabled: %s", enable_sim ? "true" : "false");
	printf("\n        I/O enabled: %s", enable_sim_io ? "true" : "false");

	if(enable_sim_io) {
		printf("\n        sim_filename_pattern: \"%s\"", sim_filename_pattern);
	}

	printf("\n        N: %u", num_bodies);

	#ifdef PAIRWISE_SOLVER_DECIMATION
	printf("\n        PAIRWISE_SOLVER_DECIMATION: %u", (unsigned int)PAIRWISE_SOLVER_DECIMATION);
	#endif

	printf("\n    Visualization:");
	printf("\n        enabled: %s", enable_vis ? "true" : "false");
	printf("\n        I/O enabled: %s", enable_vis_io ? "true" : "false");

	if(enable_vis_io) {
		printf("\n        vis_filename_pattern: \"%s\"", vis_filename);
	}

	printf("\n        dimensions: %ux%u", visualization_dimensions[0], visualization_dimensions[1]);
	printf("\n\n");

	fflush(stdout);

	if(enable_vis) {
		if( STARFLOOD_SUCCESS != visualization_init(&vis, visualization_dimensions[0], visualization_dimensions[1]) ) {
			fprintf(stderr, "fatal error: %s failed.\n", "visualization_init()");

			return EXIT_FAILURE;
		}
	}

	if(enable_sim || enable_sim_io) {
		if( STARFLOOD_SUCCESS != simulation_init(&sim, num_bodies) ) {
			fprintf(stderr, "fatal error: %s failed.\n", "simulation_init()");

			if(enable_vis) {
				visualization_free(&vis);
			}

			return EXIT_FAILURE;
		}
	}

	simulation_read(&sim, "./step_0500.raw");

	#ifdef ENABLE_SIM
	for(unsigned int step_num = 0u; step_num <= num_timesteps; step_num++) {
	#else
	for(unsigned int step_num = 0u; step_num <= num_timesteps; step_num += (unsigned int)OUTPUT_INTERVAL) {
	#endif
		printf("Step #%u (run is %6.03f%% complete)\n", step_num, 100.0 * ( (double)step_num / (double)num_timesteps ) );
		printf("t = %.06f\n", (double)TIMESTEP_SIZE * (double)step_num);
		//fflush(stdout);

		sim.step_number = step_num;

		// check stop condition (a file named "stop" in OUTPUT_DIR)
		{
			FILE* stopfile = fopen(OUTPUT_DIR "/" "stop", "r");

			if(NULL != (void*)stopfile) {
				fprintf(stderr, "Stopfile found! Stopping run...\n");

				if( 0 != fclose(stopfile) ) {
					fprintf(stderr, "error: fclose(stopfile) ");
					perror("failed");
				}

				if( 0 != remove(OUTPUT_DIR "/stop") ) {
					fprintf(stderr, "error: remove(\"%s\") ", OUTPUT_DIR "/stop");
					perror("failed");
				}

				break;
			}
		}

		#ifdef SIM_FILENAME
		unsigned int sim_file_num = step_num;

		#ifdef OUTPUT_INTERVAL
		sim_file_num /= (unsigned int)OUTPUT_INTERVAL;
		#endif

		if( 0 >= snprintf(sim_filename, filename_max, sim_filename_pattern, sim_file_num) ) {
			fprintf(stderr, "fatal error: snprintf(sim_filename, filename_max, \"%s\", %u) ", sim_filename_pattern, step_num);
			perror("failed");
			break;
		}
		#endif

		#ifdef VIS_FILENAME
		unsigned int vis_file_num = step_num;

		#ifdef OUTPUT_INTERVAL
		vis_file_num /= (unsigned int)OUTPUT_INTERVAL;
		#endif

		if( 0 >= snprintf(vis_filename, filename_max, vis_filename_pattern, vis_file_num) ) {
			fprintf(stderr, "fatal error: snprintf(vis_filename, filename_max, \"%s\", %u) ", vis_filename_pattern, step_num);
			perror("failed");
			break;
		}
		#endif

		#ifdef OUTPUT_INTERVAL
		if( (0u == (step_num % OUTPUT_INTERVAL)) && enable_sim_io ) {
		#else
		if(enable_sim_io) {
		#endif
			if(enable_sim) {
				if( STARFLOOD_SUCCESS != simulation_save(&sim, sim_filename) ) {
					fprintf(stderr, "error: simulation_save(&sim, \"%s\") failed.\n", sim_filename);
				}
			} else {
				if( STARFLOOD_SUCCESS != simulation_read(&sim, sim_filename) ) {
					fprintf(stderr, "fatal error: simulation_read(&sim, \"%s\") failed.\n", sim_filename);
					break;
				}
			}
		}

		#ifdef OUTPUT_INTERVAL
		if( (0u == (step_num % OUTPUT_INTERVAL)) && enable_vis) {
		#else
		if(enable_vis) {
		#endif
			if( STARFLOOD_SUCCESS != visualization_draw(&vis, &sim) ) {
				fprintf(stderr, "error: visualization_draw(&vis, &sim) failed.\n");
			}

			if(enable_vis_io) {
				if( STARFLOOD_SUCCESS != visualization_save(&vis, vis_filename) ) {
					fprintf(stderr, "error: visualization_save(&vis, \"%s\") failed.\n", vis_filename);
				}
			}
		}

		if(enable_sim) {
			if( STARFLOOD_SUCCESS != simulation_step(&sim) ) {
				fprintf(stderr, "error: simulation_step(&sim) failed.\n");
			}
		}

		printf("\n");
	}

	printf("Finished! Cleaning up...\n");

	fflush(stdout);

	if(enable_sim || enable_sim_io) {
		simulation_free(&sim);
	}

	if(enable_vis) {
		visualization_free(&vis);
	}

	free(sim_filename);
	free(vis_filename);
	free(sim_filename_pattern);
	free(vis_filename_pattern);

	printf("Done. Goodbye!\n");

	return EXIT_SUCCESS;
}
