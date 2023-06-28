#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>

#include "config.h"
#include "types.h"
#include "image.h"
#include "rng.h"

#ifdef STARFLOOD_ENABLE_PROFILING
#include <chrono>
#endif

//const real G = 6.6743015e-11; // 2018 CODATA-recommended value of the gravitational constant: 6.67430(15)e-11*m^3*kg^-1*s^-2
const real G = (real)1;

// particle-particle solver, O(N^2)
void particle_particle(real* a, real* x, real* m, int N) {
	#pragma omp simd
	for(int i = 0; i < 2 * N; i++) a[i] = (real)0;

	// the problematic part
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(i == j) continue; // don't calculate forces with itself

			real dx = x[2 * j + 0] - x[2 * i + 0], dy = x[2 * j + 1] - x[2 * i + 1];

			real r2 = dx*dx+dy*dy;

			a[2 * i + 0] += (m[j] * dx) / ((real)sqrtf(r2) * (r2 + (real)0.01));
			a[2 * i + 1] += (m[j] * dy) / ((real)sqrtf(r2) * (r2 + (real)0.01));
		}
	}

	#pragma omp simd
	for(int i = 0; i < 2 * N; i++) a[i] *= G * m[i / 2];
}

int main(int argc, char** argv) {
	real* x; // position buffer
	real* v; // velocity buffer
	real* a; // acceleration buffer
	real* m; // mass buffer
	float* image; // image buffer

	int N = 1000; // number of bodies
	real dt = (real)0.1; // timestep size
	int image_size_x = 512, image_size_y = 512; // image size

	if(argc > 1) {
		printf("%s: invalid option \'%s\'\n", argv[0], argv[1]);

		return EXIT_FAILURE;
	}

	printf(
	"+---------------+\n"
	"|   Starflood   |\n"
	"+---------------+\n"
	"\n"
	"Parameters:\n"
	"{\n"
	"  N: %d\n"
	"  dt: %f\n"
	"  Image Size: %dx%d\n"
	"}\n\n", N, dt, image_size_x, image_size_y);

	fflush(stdout);

	size_t x_size = (size_t)N * (size_t)2 * sizeof(real);
	size_t v_size = (size_t)N * (size_t)2 * sizeof(real);
	size_t a_size = (size_t)N * (size_t)2 * sizeof(real);
	size_t m_size = (size_t)N * (size_t)1 * sizeof(real);
	size_t t_size = x_size + v_size + a_size + m_size;
	size_t image_size = (size_t)image_size_x * (size_t)image_size_y * (size_t)3 * sizeof(float); // RGB32F

	// allocate memory
	{
		printf("Allocating %zu B of memory for simulation...", t_size);

		fflush(stdout);

		x = (real*)malloc(x_size);

		if(x == NULL) {
			perror(" malloc() returned NULL when allocating x! Error");

			return EXIT_FAILURE;
		}

		v = (real*)malloc(v_size);

		if(v == NULL) {
			perror(" malloc() returned NULL when allocating v! Error");

			free(x);

			return EXIT_FAILURE;
		}

		a = (real*)malloc(a_size);

		if(a == NULL) {
			perror(" malloc() returned NULL when allocating a! Error");

			free(v);
			free(x);

			return EXIT_FAILURE;
		}

		m = (real*)malloc(m_size);

		if(m == NULL) {
			perror(" malloc() returned NULL when allocating m! Error");

			free(a);
			free(v);
			free(x);

			return EXIT_FAILURE;
		}

		printf(" Done!\n");

		printf("Allocating %zu B of memory for rendering...", image_size);

		fflush(stdout);

		image = (float*)malloc(image_size);

		printf(" Done!\n");
	}

	// initialize simulation
	{
		rng_state_t ic_rng_state; // initial condition rng state

		ic_rng_state.ns = (uint32_t)74u;
		update_rng(&ic_rng_state);

		// initialize position
		//#pragma omp simd
		for(int i = 0; i < N; i++) {
			x[2 * i + 0] = (real)(20.0f * (urand(&ic_rng_state) - 0.5f));
			x[2 * i + 1] = (real)(20.0f * (urand(&ic_rng_state) - 0.5f));
		}

		ic_rng_state.ns = (uint32_t)67u;
		update_rng(&ic_rng_state);

		// initialize velocity
		//#pragma omp simd
		for(int i = 0; i < N; i++) {
			//v[2 * i + 0] = (real)0;
			//v[2 * i + 1] = (real)0;
			v[2 * i + 0] = (real)(0.02f * (urand(&ic_rng_state) - 0.5f));
			v[2 * i + 1] = (real)(0.02f * (urand(&ic_rng_state) - 0.5f));
		}

		ic_rng_state.ns = (uint32_t)14u;
		update_rng(&ic_rng_state);

		// initialize mass
		{
			real m0 = 10.0 / (double)N;

			#pragma omp simd
			for(int i = 0; i < N; i++) m[i] = m0;
		}

		particle_particle(a, x, m, N); // update acceleration
	}

	#ifdef STARFLOOD_ENABLE_PROFILING
	std::chrono::time_point<std::chrono::high_resolution_clock> t0, t1;

	FILE* diagfile = fopen("./log.csv", "w");

	fprintf(diagfile, "\"n\",\"Render Time (ms)\",\"First 1/2 Kick (ms)\",\"Drift (ms)\",\"Update Acceleration (ms)\",\"Second 1/2 Kick (ms)\"\n");

	fflush(diagfile);
	#endif

	// run and render simulation
	for(int n = 0; n < 100; n++) {
		#ifdef STARFLOOD_ENABLE_PROFILING
		fprintf(diagfile, "%d,", n);
		#endif

		// render simulation
		//if(n % FRAME_INTERVAL == 0) {
		if(true) {
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = std::chrono::high_resolution_clock::now();
			#endif

			#pragma omp simd
			for(int i = 0; i < (3 * image_size_x * image_size_y); i++) image[i] = 0.0f;

			for(int i = 0; i < N; i++) {
				// vec2 uv = (fragCoord - 0.5 * resolution) / resolution.y
				real uv[2] = {x[2 * i + 0], x[2 * i + 1]};

				uv[0] *= (real)0.025; // 0.05
				uv[1] *= (real)0.025;

				int coord[2] = {
				(int)((real)image_size_y * uv[0] + (real)0.5 * (real)image_size_x),
				(int)((real)image_size_y * uv[1] + (real)0.5 * (real)image_size_y)};

				if(0 <= coord[0] && coord[0] < image_size_x && 0 <= coord[1] && coord[1] < image_size_y) {
					image[3 * (image_size_x * coord[1] + coord[0]) + 0] += 0.5f;
					image[3 * (image_size_x * coord[1] + coord[0]) + 1] += 0.5f;
					image[3 * (image_size_x * coord[1] + coord[0]) + 2] += 0.5f;
				}
			}

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = std::chrono::high_resolution_clock::now();

			fprintf(diagfile, "%f,", (double)((std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0)).count()) / 1000000.0);
			#endif

			char filename[sizeof "./out1/img_0000.hdr"];

			sprintf(filename, "./out/img_%04d.hdr", n);

			write_image_hdr(filename, image_size_x, image_size_y, image);
		}
		#ifdef STARFLOOD_ENABLE_PROFILING
		else {
			fprintf(diagfile, "0.0,");
		}
		#endif

		#ifdef STARFLOOD_ENABLE_PROFILING
		t0 = std::chrono::high_resolution_clock::now();
		#endif
		// VERY FAST
		#pragma omp simd
		for(int i = 0; i < 2 * N; i++) v[i] += (real)0.5 * dt * a[i]; // (1/2) kick
		#ifdef STARFLOOD_ENABLE_PROFILING
		t1 = std::chrono::high_resolution_clock::now();
		fprintf(diagfile, "%f,", (double)((std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0)).count()) / 1000000.0);
		#endif

		#ifdef STARFLOOD_ENABLE_PROFILING
		t0 = std::chrono::high_resolution_clock::now();
		#endif
		// VERY FAST
		#pragma omp simd
		for(int i = 0; i < 2 * N; i++) x[i] += dt * v[i]; // drift
		#ifdef STARFLOOD_ENABLE_PROFILING
		t1 = std::chrono::high_resolution_clock::now();
		fprintf(diagfile, "%f,", (double)((std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0)).count()) / 1000000.0);
		#endif

		#ifdef STARFLOOD_ENABLE_PROFILING
		t0 = std::chrono::high_resolution_clock::now();
		#endif
		// VERY SLOW
		particle_particle(a, x, m, N); // update acceleration
		#ifdef STARFLOOD_ENABLE_PROFILING
		t1 = std::chrono::high_resolution_clock::now();
		fprintf(diagfile, "%f,", (double)((std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0)).count()) / 1000000.0);
		#endif

		#ifdef STARFLOOD_ENABLE_PROFILING
		t0 = std::chrono::high_resolution_clock::now();
		#endif
		// VERY FAST
		#pragma omp simd
		for(int i = 0; i < 2 * N; i++) v[i] += (real)0.5 * dt * a[i]; // (1/2) kick
		#ifdef STARFLOOD_ENABLE_PROFILING
		t1 = std::chrono::high_resolution_clock::now();
		fprintf(diagfile, "%f\n", (double)((std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0)).count()) / 1000000.0);
		#endif
	}

	#ifdef STARFLOOD_ENABLE_PROFILING
	fclose(diagfile);
	#endif

	// clean up
	{
		printf("Freeing %zu B of memory...", image_size + t_size);

		fflush(stdout);

		free(image);
		free(m);
		free(a);
		free(v);
		free(x);

		printf(" Done!\n");
	}

	printf("Finished!\n");

	return EXIT_SUCCESS;
}
