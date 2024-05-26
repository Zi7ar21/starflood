#include <simulation/simulation.hpp>

#include <cstdio>
#include <cstdlib>

void sim::clean() {
	const real dt = TIMESTEP; // Simulation Timestep

	// Clean up
	{
		printf("\r\x1b[2K"); // return to beginning and clear line (ANSI escape sequence used for clearing line)

		printf("[%.09f] Run complete!\n", omp_get_wtime()-start_time);

		printf("[%.09f] Cleaning up.\n", omp_get_wtime()-start_time);

		fflush(stdout);

		#ifdef STARFLOOD_ENABLE_PROFILING
		fclose(diagfile);
		#endif

		#ifdef STARFLOOD_ENABLE_STATISTICS
		fclose(statistics_file);
		#endif

		printf("[%.09f] Attempting to free %zu bytes (%.03f MiB) used during the run.\n", omp_get_wtime()-start_time, mem_size, 0.00000095367431640625*(double)mem_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		free(mem);

		t1 = omp_get_wtime();

		printf("[%.09f] %zu bytes (%.03f MiB) were freed (took %.9f seconds).\n", omp_get_wtime()-start_time, mem_size, 0.00000095367431640625*(double)mem_size, t1 - t0);

		printf("[%.09f] Attempting to free %zu bytes (%.03f MiB) used during rendering.\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		free(image);

		t1 = omp_get_wtime();

		printf("[%.09f] %zu bytes (%.03f MiB) were freed (took %.9f seconds).\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size, t1 - t0);

		printf("[%.09f] Execution complete. Thank you for using Starflood!\n", omp_get_wtime()-start_time);

		fflush(stdout);
	}
}
