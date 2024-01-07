#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstdint>
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <common.h>
#include <rng.hpp>
#include <barnes-hut.hpp>

void drawLine(float* image, int w, int h, float r, float g, float b, float a, int x0, int y0, int x1, int y1) {
	// Bresenham's Line Algorithm

	int dx = abs(x1 - x0);
	int sx = x0 < x1 ? 1 : -1;
	int dy = -abs(y1 - y0);
	int sy = y0 < y1 ? 1 : -1;
	int error = dx + dy;

	while(true) {
		if((0 <= x0) && (x0 < w) && (0 <= y0) && (y0 < h)) {
			int index = (w*y0)+x0;
			image[4*index+0] += r;
			image[4*index+1] += g;
			image[4*index+2] += b;
			image[4*index+3] += a;
		}

		if(x0 == x1 && y0 == y1) break;

		int e2 = 2 * error;

		if(e2 >= dy) {
			if(x0 == x1) break;
			error = error + dy;
			x0 = x0 + sx;
		}

		if(e2 <= dx) {
			if(y0 == y1) break;
			error = error + dx;
			y0 = y0 + sy;
		}
	}
}

void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1) {
	drawLine(image, w, h, r, g, b, a,
	(int)((real)h * x0 + (real)0.5 * (real)w),
	(int)((real)h * y0 + (real)0.5 * (real)h),
	(int)((real)h * x1 + (real)0.5 * (real)w),
	(int)((real)h * y1 + (real)0.5 * (real)h));
}





int main(int argc, char** argv) {
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

	printf("Starflood %d.%d.%d\n", STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH); // print version information

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

	int num_devices = omp_get_num_devices();

	printf("OpenMP found %d devices.\n", num_devices);

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

	// timestamps
	double t0, t1;

	// call omp_get_wtime() so it is cached or something idk
	t0 = omp_get_wtime();
	t1 = omp_get_wtime();

	size_t w = 640; // image width
	size_t h = 360; // image height
	size_t N = 16384; // number of bodies in the simulation

	if(N < 1) {
		printf("Error: N was set to %zu, it must be greater than 0!\n", N);

		return EXIT_FAILURE;
	}

	const real dt = TIMESTEP; // Simulation Timestep

	float* image; // Image Buffer, RGBA32F
	real* sim; // Simulation Buffer
	real* mas; // Mass Buffer
	real* pos; // Position Buffer
	real* vel; // Velocity Buffer
	real* acc; // Acceleration Buffer
	int* ids; // id buffer, contains ids of node paticles are in

	stbi_flip_vertically_on_write(1);
	size_t image_size = sizeof(float)*(4*w*h); // Image Buffer Size, RGBA32F
	size_t mas_size = (N * 1) * sizeof(real); // mass buffer size
	size_t pos_size = (N * 2) * sizeof(real); // position buffer size
	size_t vel_size = (N * 2) * sizeof(real); // velocity buffer size
	size_t acc_size = (N * 2) * sizeof(real); // acceleration buffer size
	size_t sim_size = mas_size+pos_size+vel_size+acc_size; // simulation total buffer size
	size_t ids_size = (N * 1) * sizeof(int);

	std::vector<Node> tree;

	// Allocate Memory
	{
		image = (float*)malloc(image_size); // Allocate image buffer

		if(image == NULL) {
			printf("Error: malloc() returned a null pointer!\n");

			fflush(stdout);

			perror("Error");

			return EXIT_FAILURE;
		}

		sim = (real*)malloc(sim_size); // Allocate simulation memory

		if(sim == NULL) {
			printf("Error: malloc() returned a null pointer!\n");

			fflush(stdout);

			perror("Error");

			free(image);

			return EXIT_FAILURE;
		}

		ids = (int*)malloc(ids_size);

		mas = sim; // mass buffer
		pos = mas+(mas_size/sizeof(real)); // position buffer
		vel = pos+(pos_size/sizeof(real)); // velocity buffer
		acc = vel+(vel_size/sizeof(real)); // acceleration buffer
	}

	// Initialize Simulation
	{
		for(size_t i = 0; i < (N * 1); i++) mas[i] = (real)(((long double)10)/((long double)N));
		for(size_t i = 0; i < (N * 2); i++) acc[i] = (real)0;
		for(size_t i = 0; i < (N * 2); i++) vel[i] = (real)0;

		// randomly selected positions and velocities
		for(size_t i = 0; i < N; i++) {
			uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			// normal distribution
			//pos[2*i+0] = 0.1*sqrt(-2.0*log(z0))*cos(TAU*z2);
			//pos[2*i+1] = 0.1*sqrt(-2.0*log(z0))*sin(TAU*z2);
			//pos[2*i+2] = 0.1*sqrt(-2.0*log(z1))*cos(TAU*z3);
			pos[2*i+0] = 0.2*sqrt(z1)*cos(TAU*z0);
			pos[2*i+1] = 0.2*sqrt(z1)*sin(TAU*z0);
		}

		for(size_t i = 0; i < N; i++) {
			vel[2*i+0] = -.3*pos[2*i+1];
			vel[2*i+1] =  .3*pos[2*i+0];
		}

		// calculate initial gravitational accelerations
		BarnesHut(tree, image, w, h, ids, mas, pos, acc, N, 0);
	}

	#ifdef STARFLOOD_ENABLE_PROFILING
	FILE* diagfile = fopen("./log.csv", "w");

	fprintf(diagfile, "\"n\",\"Clear Image (ms)\",\"Render Bodies (ms)\",\"Render Quadtree (ms)\",\"Write Image (ms)\",\"First Kick (ms)\",\"Drift (ms)\",\"Update Acceleration (ms)\",\"Second Kick (ms)\"\n");

	fflush(diagfile);
	#endif

	const int num_steps = FRAMES;

	for(int step_num = 0; step_num < num_steps; step_num++) {
		printf("\rRunning simulation, %d/%d (%.2f%%) completed...", step_num, num_steps, 100.0*((double)step_num/(double)num_steps));

		fflush(stdout);

		#ifdef STARFLOOD_ENABLE_PROFILING
		fprintf(diagfile, "%d,", step_num);
		#endif

		// render and write the image
		{
			stbi_flip_vertically_on_write(1);

			// clear the image buffer
			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0) {
					// VERY FAST
					//#pragma omp simd
					for(size_t i = 0; i < (4*w*h); i++) image[i] = (float)0;
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			// render the image
			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0){
					for(size_t i = 0; i < N; i++) {
						// vec2 uv = (fragCoord - 0.5 * resolution) / resolution.y
						real t = 0.04*(real)step_num;
						//t = 0.0;
						//real uv[2] = {r[2 * i + 0]*cos(t)+r[2 * i + 2]*sin(t), r[2 * i + 1]};
						real uv[2] = {pos[2*i+0], pos[2*i+1]};

						uv[0] *= (real)1.0; // 0.05
						uv[1] *= (real)1.0;

						int coord[2] = {
						(int)((real)h * uv[0] + (real)0.5 * (real)w),
						(int)((real)h * uv[1] + (real)0.5 * (real)h)};

						if((0 <= coord[0]) && ((size_t)coord[0] < w) && (0 <= coord[1]) && ((size_t)coord[1] < h)) {
							image[4*(w*coord[1]+coord[0])+0] += 0.05f;
							image[4*(w*coord[1]+coord[0])+1] += 0.05f;
							image[4*(w*coord[1]+coord[0])+2] += 0.05f;
							image[4*(w*coord[1]+coord[0])+3] += 0.0f;
						}
					}

					for(size_t i = 0; i < (w * h); i++) {
						image[4*i+3] = (real)1;
					}
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0) {
					/*
					for(size_t i = 0; i < tree.size(); i++) {
						real x_min = tree[i].x_min, y_min = tree[i].y_min, x_max = tree[i].x_max, y_max = tree[i].y_max;

						if((x_max-x_min) < (real)0.1) continue;

						drawLineUV(image, w, h, 0., .01, 0., .1, x_min, y_min, x_max, y_min);
						drawLineUV(image, w, h, 0., .01, 0., .1, x_min, y_min, x_min, y_max);
						drawLineUV(image, w, h, 0., .01, 0., .1, x_min, y_max, x_max, y_max);
						drawLineUV(image, w, h, 0., .01, 0., .1, x_max, y_min, x_max, y_max);
					}
					*/
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			if((step_num % FRAME_INTERVAL) == 0){
				#ifdef STARFLOOD_RENDER_INTERACTS
				BarnesHut(tree, image, w, h, ids, mas, pos, acc, N, step_num);
				#endif

				// test, test, is this thing working?
				drawLine(image, w, h, 1., 0., 0., 1., 15+0, 15-1, 15+4, 15-1); drawLine(image, w, h, 1., 0., 0., 1., 15+0, 15+1, 15+4, 15+1);
				drawLine(image, w, h, 1., 1., 0., 1., 15-1, 15+0, 15-1, 15+4); drawLine(image, w, h, 1., 1., 0., 1., 15+1, 15+0, 15+1, 15+4);
				drawLine(image, w, h, 0., 1., 0., 1., 15+0, 15-1, 15-4, 15-1); drawLine(image, w, h, 0., 1., 0., 1., 15+0, 15+1, 15-4, 15+1);
				drawLine(image, w, h, 0., 0., 1., 1., 15-1, 15+0, 15-1, 15-4); drawLine(image, w, h, 0., 0., 1., 1., 15+1, 15+0, 15+1, 15-4);

				for(size_t i = 0; i < (w * h); i++) {
					image[4*i+0] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
					image[4*i+1] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
					image[4*i+2] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
					image[4*i+3]  = (float)1;
				}
			}

			// write the image
			{
				char file_name[FILENAME_MAX];

				snprintf(file_name, FILENAME_MAX, "out/step_%04d.hdr", step_num/FRAME_INTERVAL);

				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0){
					#ifndef STARFLOOD_DISABLE_IMAGE_WRITE
					if(stbi_write_hdr(file_name, w, h, 4, image) == 0) {
						printf("Error: stbi_write_hdr(%s, %zu, %zu, %d, %p) was not successful!\n", file_name, w, h, 4, image);
					}
					#endif
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}
		}

		// kick
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			for(size_t i = 0; i < (N * 2); i++) vel[i] += (real)0.5 * dt * acc[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
			#endif
		}

		// drift
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			for(size_t i = 0; i < (N * 2); i++) pos[i] += dt * vel[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
			#endif
		}

		// update acceleration
		{
			for(size_t i = 0; i < (N * 2); i++) acc[i] = (real)0;

			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			BarnesHut(tree, image, w, h, ids, mas, pos, acc, N, step_num);

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
			#endif
		}

		//for(int i = 0; i < N; i++) {
		//	a[2*i+0] += -1.0*r[2*i+2];
		//	a[2*i+2] +=  1.0*r[2*i+0];
		//}

		// kick
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			//#pragma omp target teams distribute parallel for map(tofrom:vel[0:2*N]) map(tofrom:acc[0:2*N])
			for(size_t i = 0; i < (N * 2); i++) vel[i] += (real)0.5 * dt * acc[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f\n", 1000.0*(t1-t0));
			#endif
		}

		/*
		for(int i = 0; i < N; i++) {
			if(abs(pos[2*i+0] > 1.) || abs(pos[2*i+0]) > 1.) {
				// Generate Initial Conditions
				uint32_t ns = (uint32_t)1u;

				// randomly selected positions and velocities
				for(int i = 0; i < N; i++) {
					ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

					float z0 = urand(&ns);
					float z1 = urand(&ns);
					float z2 = urand(&ns);
					float z3 = urand(&ns);

					// normal distribution
					pos[2*i+0] = 0.1*sqrt(-2.0*log(z0))*cos(TAU*z2);
					pos[2*i+1] = 0.1*sqrt(-2.0*log(z0))*sin(TAU*z2);
					//pos[2*i+2] = 0.1*sqrt(-2.0*log(z1))*cos(TAU*z3);
				}
			}
		}
		*/
	}

	printf("\rRunning simulation, %d/%d (100.00%%) completed... Done!\n\n", num_steps, num_steps);

	fflush(stdout);

	#ifdef STARFLOOD_ENABLE_PROFILING
	fclose(diagfile);
	#endif

	free(image);
	free(sim);

	return EXIT_SUCCESS;
}
