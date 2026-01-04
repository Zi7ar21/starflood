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
#include "initcond.h"
#include "simulation.h"
#include "solver.h"
#include "timing.h"
#include "types.h"
#include "visualization.h"

#if (0 >= SIM_FILE_FORMAT_READ)
#define SIM_FILE_EXTENSION_READ ".raw"
#endif

#if (1 <= SIM_FILE_FORMAT_READ)
#define SIM_FILE_EXTENSION_READ ".ply"
#endif

#if (0 >= SIM_FILE_FORMAT_SAVE)
#define SIM_FILE_EXTENSION_SAVE ".raw"
#endif

#if (1 <= SIM_FILE_FORMAT_SAVE)
#define SIM_FILE_EXTENSION_SAVE ".ply"
#endif

#if (0 == VIS_FILE_FORMAT)
#define VIS_FILE_EXTENSION ".pfm"
#endif

#if (1 == VIS_FILE_FORMAT)
#define VIS_FILE_EXTENSION ".ppm"
#endif

#if (2 == VIS_FILE_FORMAT)
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
	unsigned int visualization_dimensions[2] = {1080u, 1080u}; //  Full High-Definition
	//unsigned int visualization_dimensions[2] = {2160u, 2160u}; // Ultra High-Definition

	// (16:9) Aspect Ratio
	//unsigned int visualization_dimensions[2] = { 640u,  360u}; //   Standard Definition
	//unsigned int visualization_dimensions[2] = {1280u,  720u}; //       High-Definition
	//unsigned int visualization_dimensions[2] = {1920u, 1080u}; //  Full High-Definition
	//unsigned int visualization_dimensions[2] = {3840u, 2160u}; // Ultra High-Definition

	// Cinematic Aspect Ratio
	//unsigned int visualization_dimensions[2] = {2048u, 1080u}; // DCI 2K
	//unsigned int visualization_dimensions[2] = {4096u, 2160u}; // DCI 4K

	sim_t sim;
	vis_t vis;

	#ifdef STARFLOOD_FILENAME_MAX
	const size_t filename_max      = (size_t)STARFLOOD_FILENAME_MAX;
	#else
	const size_t filename_max      = (size_t)FILENAME_MAX;
	#endif
	const size_t filename_mem_size = sizeof(char) * filename_max * (size_t)6u;

	char* filename_mem = (char*)malloc(filename_mem_size);

	if(NULL == (void*)filename_mem) {
		fprintf(stderr, "fatal error: filename_mem is (char*)NULL after malloc(%zu", filename_mem_size);
		perror(")");
		return EXIT_FAILURE;
	}

	memset(filename_mem, 0, filename_mem_size);

	char* sim_read_file_pattern = &filename_mem[(size_t)0u * filename_max];
	char* sim_save_file_pattern = &filename_mem[(size_t)1u * filename_max];
	char* vis_save_file_pattern = &filename_mem[(size_t)2u * filename_max];
	char* sim_read_file_name    = &filename_mem[(size_t)3u * filename_max];
	char* sim_save_file_name    = &filename_mem[(size_t)4u * filename_max];
	char* vis_save_file_name    = &filename_mem[(size_t)5u * filename_max];

	#ifdef SIM_FILENAME
	strcpy(sim_read_file_pattern, OUTPUT_DIR "/" SIM_FILENAME SIM_FILE_EXTENSION_READ);
	strcpy(sim_save_file_pattern, OUTPUT_DIR "/" SIM_FILENAME SIM_FILE_EXTENSION_SAVE);
	#else
	strcpy(sim_read_file_pattern, OUTPUT_DIR "/sim/step_%04u" SIM_FILE_EXTENSION_READ);
	strcpy(sim_save_file_pattern, OUTPUT_DIR "/sim/step_%04u" SIM_FILE_EXTENSION_SAVE);
	#endif

	#ifdef VIS_FILENAME
	strcpy(vis_save_file_pattern, OUTPUT_DIR "/" VIS_FILENAME VIS_FILE_EXTENSION);
	#else
	strcpy(vis_save_file_pattern, OUTPUT_DIR "/vis/step_%04u" VIS_FILE_EXTENSION);
	#endif

	if( 0 >= snprintf(sim_read_file_name, filename_max, sim_read_file_pattern, 0u) ) {
		fprintf(stderr, "fatal error: snprintf(%s, filename_max, \"%s\", %u) ", "sim_read_file_name", sim_read_file_pattern, 0u);
		perror("failed");
		free(filename_mem);
		return EXIT_FAILURE;
	}

	if( 0 >= snprintf(sim_save_file_name, filename_max, sim_save_file_pattern, 0u) ) {
		fprintf(stderr, "fatal error: snprintf(%s, filename_max, \"%s\", %u) ", "sim_save_file_name", sim_save_file_pattern, 0u);
		perror("failed");
		free(filename_mem);
		return EXIT_FAILURE;
	}

	if( 0 >= snprintf(vis_save_file_name, filename_max, vis_save_file_pattern, 0u) ) {
		fprintf(stderr, "fatal error: snprintf(%s, filename_max, \"%s\", %u) ", "vis_save_file_name", vis_save_file_pattern, 0u);
		perror("failed");
		free(filename_mem);
		return EXIT_FAILURE;
	}

	// Parse arguments
	for(int i = 0; i < argc; i++) {
		if(0 >= i) {
			continue;
		}

		if(NULL == (void*)argv[i]) {
			fprintf(stderr, "%s: error: argv[%d] is NULL!\n", argv[0], i);

			if(1 >= i) {
				break;
			}

			free(filename_mem);
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
				free(filename_mem);
				return EXIT_SUCCESS;
			}

			continue;
		}

		if( 0 == strcmp("-V", argv[i]) || 0 == strcmp("--version", argv[i]) ) {
			printf(
				"Starflood version %d.%d.%d\n"
				"\n"
				"Starflood repository: <%s>\n",
				STARFLOOD_VERSION_MAJOR,
				STARFLOOD_VERSION_MINOR,
				STARFLOOD_VERSION_PATCH,
				STARFLOOD_REPOSITORY
			);

			if( (argc - 1) == i ) {
				free(filename_mem);
				return EXIT_SUCCESS;
			}

			continue;
		}

		fprintf(stderr, "%s: error: unrecognized option \'%s\'\nTry \'%s --help\' for more information.\n", argv[0], argv[1], argv[0]);

		free(filename_mem);
		return EXIT_FAILURE;
	}

	printf("=== Starflood ===\n");

	printf("  Version %d.%d.%d\n", STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH);

	fflush(stdout);

	printf("\nConfiguration:");
	printf("\n    OpenMP Support:");
	printf("\n        Enabled: ");
	#ifdef _OPENMP
	printf("true");
	printf("\n        _OPENMP: %d", _OPENMP);
	printf("\n        max_threads: %d", omp_get_max_threads());
	printf("\n        num_devices: %d", omp_get_num_devices());
	printf("\n        wtick: %.03f nsec", 1.0e9*omp_get_wtick());
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
			struct timespec ts0, ts1;

			printf("\n            res: ");

			if( 0 != clock_getres(STARFLOOD_POSIX_CLOCKID, &ts0) ) {
				printf("N/A\n");
			} else {
				printf("%jd nsec", (intmax_t)1000000000l*(intmax_t)ts0.tv_sec+(intmax_t)ts0.tv_sec);
			}

			printf("\n            time: ");

			if( 0 != clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts0) ) {
				printf("N/A\n");
			} else {
				printf("%jd.%09jd sec", (intmax_t)ts0.tv_sec, (intmax_t)ts0.tv_nsec);
			}

			printf("\n            minimum difference (100 samples): ");

			intmax_t clock_minimum_difference = 0;

			// first sample is discarded
			for(int i = 0; i <= 100; i++) {
				// call clock_gettime multiple times for hot caching
				// technically this entire thing is kind of pointless
				// because of processor frequency scaling
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts0);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts1);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts0);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts1);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts0);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts1);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts0);
				clock_gettime(STARFLOOD_POSIX_CLOCKID, &ts1);

				clock_minimum_difference += (0 < i) ? (intmax_t)1000000000l*(intmax_t)(ts1.tv_sec-ts0.tv_sec)+((intmax_t)ts1.tv_nsec-(intmax_t)ts0.tv_nsec) : (intmax_t)0;
			}

			printf("%.03f nsec", 0.01 * (double)clock_minimum_difference);
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
		printf("\n        sim_read_file_pattern: \"%s\"", sim_read_file_pattern);
		printf("\n        sim_save_file_pattern: \"%s\"", sim_save_file_pattern);
	}

	printf("\n        N: %u", num_bodies);

	#ifdef PAIRWISE_SOLVER_DECIMATION
	printf("\n        PAIRWISE_SOLVER_DECIMATION: %u", (unsigned int)PAIRWISE_SOLVER_DECIMATION);
	#endif

	printf("\n    Visualization:");
	printf("\n        enabled: %s", enable_vis ? "true" : "false");
	printf("\n        I/O enabled: %s", enable_vis_io ? "true" : "false");

	if(enable_vis_io) {
		printf("\n        vis_save_file_pattern: \"%s\"", vis_save_file_pattern);
	}

	printf("\n        dimensions: %ux%u", visualization_dimensions[0], visualization_dimensions[1]);
	printf("\n\n");

	fflush(stdout);

	// initialize the visualization
	if(enable_vis) {
		if( STARFLOOD_SUCCESS != visualization_init(&vis, visualization_dimensions[0], visualization_dimensions[1]) ) {
			fprintf(stderr, "fatal error: visualization_init(&vis, %u, %u) failed.\n", visualization_dimensions[0], visualization_dimensions[1]);
			free(filename_mem);
			return EXIT_FAILURE;
		}
	}

	// initialize the simulation
	if( STARFLOOD_SUCCESS != sim_init(&sim, num_bodies) ) {
		fprintf(stderr, "fatal error: sim_init(&sim, %u) failed.\n", num_bodies);

		if(enable_vis) {
			visualization_free(&vis);
		}

		free(filename_mem);
		return EXIT_FAILURE;
	}

	// simulation initial conditions
	#ifdef INIT_COND_FILE
	if( STARFLOOD_SUCCESS != sim_read_raw(&sim, INIT_COND_FILE) ) {
		fprintf(stderr, "fatal error: sim_read_raw(&sim, \"%s\") failed.\n", INIT_COND_FILE);

		if(enable_vis) {
			visualization_free(&vis);
		}

		sim_free(&sim);
		free(filename_mem);
		return EXIT_FAILURE;
	}
	#else
	if(enable_sim) {
		if( STARFLOOD_SUCCESS != initcond_generate(&sim) ) {
			fprintf(stderr, "fatal error: initcond_generate() failed.\n");

			if(enable_vis) {
				visualization_free(&vis);
			}

			sim_free(&sim);
			free(filename_mem);
			return EXIT_FAILURE;
		}

		// An initial solver step is needed for the first leapfrog "kick" in the "kick-drift-kick" form, since acceleration is only updated after "drift"
		if( STARFLOOD_SUCCESS != solver_run(&sim, 0u) ) {
			fprintf(stderr, "fatal error: solver_run() failed.\n");

			if(enable_vis) {
				visualization_free(&vis);
			}

			sim_free(&sim);
			free(filename_mem);
			return EXIT_FAILURE;
		}
	}
	#endif

	unsigned int step_inc = 1u;

	#ifdef OUTPUT_INTERVAL
	if(!enable_sim) {
		step_inc = OUTPUT_INTERVAL;
	}
	#endif

	for(unsigned int step_num = 0u; step_num <= num_timesteps; step_num += step_inc) {
		printf("Step #%u (run is %6.03f%% complete)\n", step_num, 100.0 * (0 < num_timesteps ? (double)step_num / (double)num_timesteps : 1.0) );
		printf("t = %.06f\n", (double)TIMESTEP_SIZE * (double)step_num);

		// check stop condition (a file named "stop" in OUTPUT_DIR)
		{
			FILE* stopfile = fopen(OUTPUT_DIR "/" "stop", "r");

			if(NULL != (void*)stopfile) {
				printf("\nStopfile found! Stopping run...\n");
				fflush(stdout);

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

		fflush(stdout);

		sim.step_number = step_num;

		#ifdef SIM_FILENAME
		unsigned int sim_file_num = step_num;

		#ifdef OUTPUT_INTERVAL
		sim_file_num /= (unsigned int)OUTPUT_INTERVAL;
		#endif

		if( 0 >= snprintf(sim_read_file_name, filename_max, sim_read_file_pattern, sim_file_num) ) {
			fprintf(stderr, "fatal error: snprintf(%s, filename_max, \"%s\", %u) ", "sim_read_file_name", sim_read_file_pattern, step_num);
			perror("failed");
			break;
		}

		if( 0 >= snprintf(sim_save_file_name, filename_max, sim_save_file_pattern, sim_file_num) ) {
			fprintf(stderr, "fatal error: snprintf(%s, filename_max, \"%s\", %u) ", "sim_save_file_name", sim_save_file_pattern, step_num);
			perror("failed");
			break;
		}
		#endif

		#ifdef VIS_FILENAME
		unsigned int vis_file_num = step_num;

		#ifdef OUTPUT_INTERVAL
		vis_file_num /= (unsigned int)OUTPUT_INTERVAL;
		#endif

		if( 0 >= snprintf(vis_save_file_name, filename_max, vis_save_file_pattern, vis_file_num) ) {
			fprintf(stderr, "fatal error: snprintf(%s, filename_max, \"%s\", %u) ", "vis_save_file_name", vis_save_file_pattern, step_num);
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
				#if (0 == SIM_FILE_FORMAT_SAVE)
				if( STARFLOOD_SUCCESS != sim_save_raw(&sim, sim_save_file_name) ) {
					fprintf(stderr, "error: sim_save_raw(&sim, \"%s\") failed.\n", sim_save_file_name);
				}
				#else
				if( STARFLOOD_SUCCESS != sim_save_ply(&sim, sim_save_file_name) ) {
					fprintf(stderr, "error: sim_save_ply(&sim, \"%s\") failed.\n", sim_save_file_name);
				}
				#endif
			} else {
				#if (0 == SIM_FILE_FORMAT_READ)
				if( STARFLOOD_SUCCESS != sim_read_raw(&sim, sim_read_file_name) ) {
					fprintf(stderr, "fatal error: sim_read_raw(&sim, \"%s\") failed.\n", sim_read_file_name);
					break;
				}
				#else
				if( STARFLOOD_SUCCESS != sim_read_ply(&sim, sim_read_file_name) ) {
					fprintf(stderr, "fatal error: sim_read_ply(&sim, \"%s\") failed.\n", sim_read_file_name);
					break;
				}
				#endif

				#if (SIM_FILE_FORMAT_READ != SIM_FILE_FORMAT_SAVE)
				#if (0 == SIM_FILE_FORMAT_SAVE)
				if( STARFLOOD_SUCCESS != sim_save_raw(&sim, sim_save_file_name) ) {
					fprintf(stderr, "error: sim_save_raw(&sim, \"%s\") failed.\n", sim_save_file_name);
				}
				#else
				if( STARFLOOD_SUCCESS != sim_save_ply(&sim, sim_save_file_name) ) {
					fprintf(stderr, "error: sim_save_ply(&sim, \"%s\") failed.\n", sim_save_file_name);
				}
				#endif
				#endif
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
				if( STARFLOOD_SUCCESS != visualization_save(&vis, vis_save_file_name) ) {
					fprintf(stderr, "error: visualization_save(&vis, \"%s\") failed.\n", vis_save_file_name);
				}
			}
		}

		if( (step_num < num_timesteps) && enable_sim) {
			if( STARFLOOD_SUCCESS != sim_step(&sim) ) {
				fprintf(stderr, "error: sim_step(&sim) failed.\n");
			}
		}

		printf("\n");
	}

	printf("Finished! Cleaning up...\n");

	fflush(stdout);

	if( STARFLOOD_SUCCESS != sim_free(&sim) ) {
		fprintf(stderr, "error: sim_free(&sim) failed.\n");
	}

	if(enable_vis) {
		if( STARFLOOD_SUCCESS != visualization_free(&vis) ) {
			fprintf(stderr, "error: visualization_free(&vis) failed.\n");
		}
	}

	free(filename_mem);

	printf("Done. Goodbye!\n");

	return EXIT_SUCCESS;
}
