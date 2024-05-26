#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstdint>
#include <omp.h>
#include <cstring>

//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <utilities/common.h>
#include <utilities/rng.hpp>
#include <gravity/barnes-hut.hpp>
#include <simulation/simulation.hpp>
#include <utilities/graphics.hpp>

int main(int argc, char** argv) {
	double start_time = omp_get_wtime();

	if(argc > 1) {
		if(!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")) {
			printf(
			"Usage: %s [OPTION]...\n"
			"  -h, --help       display this help and exit\n"
			//"      --version    output version information and exit\n"
			, argv[0]);

			return EXIT_SUCCESS;
		} else {
			printf("%s: unrecognized option '%s'\nTry '%s -h' or '%s --help' for more information.\n", argv[0], argv[1], argv[0], argv[0]);

			return EXIT_FAILURE;
		}
	}


	printf("[%.09f] Starflood %d.%d.%d\n", 0.0, STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH); // print version information

	#ifndef STARFLOOD_NO_SPLASH
	printf(
	"\n"
	"╭───┬─╴╭─╮┌─╮┌──┐  ╭─╮╭─╮┌─╮\n"
	"│   │  │ ││ ││  │  │ ││ ││ │\n"
	"╰─╮ │  ├─┤├─╯├─╴│  │ ││ ││ │\n"
	"│ │ │  │ ││ ││  │  │ ││ ││ │\n"
	"╰─╯ ╵  ╵ ╵╵ ╵╵  └──╰─╯╰─╯└─╯\n"
	"\n");
	#endif

	fflush(stdout);

	// OpenMP Offloading Stuff
	{
		printf("[%.09f] Setting up Starflood.\n", omp_get_wtime()-start_time);

		printf("[%.09f] OpenMP found %d processor(s).\n", omp_get_wtime()-start_time, omp_get_num_procs());

		int num_devices = omp_get_num_devices();

		printf("[%.09f] OpenMP found %d device(s) (for offloading).\n", omp_get_wtime()-start_time, num_devices);

		printf("[%.09f] NOTE: Starflood does not support OpenMP offloading quite yet. If you would like to see offloading in the future, check the Starflood repo for ways you can contribute!\n", omp_get_wtime()-start_time);

		printf("[%.09f] Starflood will perform calculations AND rendering using %d OpenMP thread(s).\n", omp_get_wtime()-start_time, omp_get_max_threads());

		fflush(stdout);

		/*
		#pragma omp target teams distribute parallel
		{
			printf("Hello, world!");

			//#pragma omp single
			//{
				//printf("Hello, world!");

				//double a = omp_get_wtime();

				//while(omp_get_wtime() <= (a + 5.0)) {
				//	continue;
				//}
			//}
		}
		*/
	}





	size_t w = RENDER_W; // image width
	size_t h = RENDER_H; // image height
	size_t N = NUM_BODIES; // number of bodies in the simulation

	if(N < 1) {
		printf("[%.09f] Error: N was set to %zu, it must be greater than 0! Starflood will abort now.\n", omp_get_wtime()-start_time, N);

		return EXIT_FAILURE;
	}

	stbi_flip_vertically_on_write(1);

	sim simulation(start_time);

	// call omp_get_wtime() so it is cached or something idk
	simulation.t0 = omp_get_wtime();
	simulation.t1 = omp_get_wtime();

	simulation.setup();

	simulation.init();

	simulation.run();

	simulation.clean();

	float* image; // Image Buffer, RGBA32F
	real* sim; // Simulation Buffer
	real* mas; // Mass Buffer
	real* pos; // Position Buffer
	real* vel; // Velocity Buffer
	real* acc; // Acceleration Buffer
	real* pen; // Potential Energy Buffer
	int* ids; // id buffer, contains ids of node paticles are in


	return EXIT_SUCCESS;
}
