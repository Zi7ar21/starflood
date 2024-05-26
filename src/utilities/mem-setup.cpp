#include <simulation/simulation.hpp>

#include <cstdio>
#include <cstdlib>

void sim::setup() {
	w = RENDER_W; // image width
	h = RENDER_H; // image height
	N = NUM_BODIES;

	sim_particles.mas_size = (N * 1) * sizeof(real); // mass buffer size
	sim_particles.pos_size = (N * 3) * sizeof(real); // position buffer size
	sim_particles.vel_size = (N * 3) * sizeof(real); // velocity buffer size
	sim_particles.acc_size = (N * 3) * sizeof(real); // acceleration buffer size
	sim_particles.pen_size = (N * 1) * sizeof(real); // potential energy buffer size
	sim_particles.ids_size = (N * 1) * sizeof(int);
	mem_size = sim_particles.mas_size+sim_particles.pos_size+sim_particles.vel_size+sim_particles.acc_size+sim_particles.pen_size;

	image_size = sizeof(float)*(4*w*h);

// = mas_size+pos_size+vel_size+acc_size+pen_size

	// Allocate Memory
		printf("[%.09f] Memory allocation starting.\n", omp_get_wtime()-start_time);

		printf("[%.09f] Trying to allocate %zu bytes (%.03f MiB) for rendering.\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		image = (float*)malloc(image_size); // Allocate image buffer

		t1 = omp_get_wtime();

		if(image == NULL) {
			perror("Error");

			printf("[%.09f] Error: malloc() returned a NULL pointer!\n", omp_get_wtime()-start_time);

			printf("[%.09f] %zu bytes (%.03f MiB) was allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, (size_t)0, 0.00000095367431640625*(double)0, t1 - t0);

			fflush(stdout);

			exit(EXIT_FAILURE);
		}

		printf("[%.09f] %zu bytes (%.03f MiB) was allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size, t1 - t0);

		printf("[%.09f] Trying to allocate %zu bytes (%.03f MiB) for the simulation.\n", omp_get_wtime()-start_time, mem_size, 0.00000095367431640625*(double)mem_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		mem = (real*)malloc(mem_size); // Allocate simulation memory

		t1 = omp_get_wtime();

		if(mem == NULL) {
			perror("Error");

			printf("[%.09f] Error: malloc() returned a NULL pointer!\n", omp_get_wtime()-start_time);

			printf("[%.09f] %zu bytes (%.03f MiB) was allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, (size_t)0, 0.00000095367431640625*(double)0, t1 - t0);

			fflush(stdout);

			exit(EXIT_FAILURE);
		}

		printf("[%.09f] %zu bytes (%.03f MiB) was allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, mem_size, 0.00000095367431640625*(double)mem_size, t1 - t0);

		printf("[%.09f] Trying to allocate an initial %zu bytes (%.03f MiB) for the Barnes-Hut tree.\n", omp_get_wtime()-start_time, sim_particles.ids_size, 0.00000095367431640625*(double)sim_particles.ids_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		sim_particles.ids = (int*)malloc(sim_particles.ids_size);

		t1 = omp_get_wtime();

		if(sim_particles.ids == NULL) {
			perror("Error");

			printf("[%.09f] Error: malloc() returned a NULL pointer!\n", omp_get_wtime()-start_time);

			printf("[%.09f] %zu bytes (%.03f MiB) was allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, (size_t)0, 0.00000095367431640625*(double)0, t1 - t0);

			fflush(stdout);

			exit(EXIT_FAILURE);
		}

		printf("[%.09f] %zu bytes (%.03f MiB) was allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, sim_particles.ids_size, 0.00000095367431640625*(double)sim_particles.ids_size, t1 - t0);

		fflush(stdout);

		// partition allocated memory for the simulation
		sim_particles.mas = (real*)mem; // mass buffer
		sim_particles.pos = sim_particles.mas+(sim_particles.mas_size/sizeof(real)); // position buffer
		sim_particles.vel = sim_particles.pos+(sim_particles.pos_size/sizeof(real)); // velocity buffer
		sim_particles.acc = sim_particles.vel+(sim_particles.vel_size/sizeof(real)); // acceleration buffer
		sim_particles.pen = sim_particles.acc+(sim_particles.acc_size/sizeof(real)); // potential energy buffer
}
