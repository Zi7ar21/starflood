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
#include <graphics.hpp>

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

	// profiling timestamps
	double t0, t1;

	// call omp_get_wtime() so it is cached or something idk
	t0 = omp_get_wtime();
	t1 = omp_get_wtime();

	size_t w = RENDER_W; // image width
	size_t h = RENDER_H; // image height
	size_t N = NUM_BODIES; // number of bodies in the simulation

	if(N < 1) {
		printf("[%.09f] Error: N was set to %zu, it must be greater than 0! Starflood will abort now.\n", omp_get_wtime()-start_time, N);

		return EXIT_FAILURE;
	}

	const real dt = TIMESTEP; // Simulation Timestep

	float* image; // Image Buffer, RGBA32F
	real* sim; // Simulation Buffer
	real* mas; // Mass Buffer
	real* pos; // Position Buffer
	real* vel; // Velocity Buffer
	real* acc; // Acceleration Buffer
	real* pen; // Potential Energy Buffer
	int* ids; // id buffer, contains ids of node paticles are in

	stbi_flip_vertically_on_write(1);
	size_t image_size = sizeof(float)*(4*w*h); // Image Buffer Size, RGBA32F
	size_t mas_size = (N * 1) * sizeof(real); // mass buffer size
	size_t pos_size = (N * 3) * sizeof(real); // position buffer size
	size_t vel_size = (N * 3) * sizeof(real); // velocity buffer size
	size_t acc_size = (N * 3) * sizeof(real); // acceleration buffer size
	size_t pen_size = (N * 1) * sizeof(real); // potential energy buffer size
	size_t sim_size = mas_size+pos_size+vel_size+acc_size+pen_size; // simulation total buffer size
	size_t ids_size = (N * 1) * sizeof(int);

	std::vector<Node> tree;

	// Allocate Memory
	{
		printf("[%.09f] Memory allocation starting.\n", omp_get_wtime()-start_time);

		printf("[%.09f] Trying to allocate %zu bytes (%.03f MiB) for rendering.\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		image = (float*)malloc(image_size); // Allocate image buffer

		t1 = omp_get_wtime();

		if(image == NULL) {
			perror("Error");

			printf("[%.09f] Error: malloc() returned a NULL pointer!\n", omp_get_wtime()-start_time);

			printf("[%.09f] %zu bytes (%.03f MiB) were allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, (size_t)0, 0.00000095367431640625*(double)0, t1 - t0);

			fflush(stdout);

			return EXIT_FAILURE;
		}

		printf("[%.09f] %zu bytes (%.03f MiB) were allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size, t1 - t0);

		printf("[%.09f] Trying to allocate %zu bytes (%.03f MiB) for the simulation.\n", omp_get_wtime()-start_time, sim_size, 0.00000095367431640625*(double)sim_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		sim = (real*)malloc(sim_size); // Allocate simulation memory

		t1 = omp_get_wtime();

		if(sim == NULL) {
			perror("Error");

			printf("[%.09f] Error: malloc() returned a NULL pointer!\n", omp_get_wtime()-start_time);

			printf("[%.09f] %zu bytes (%.03f MiB) were allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, (size_t)0, 0.00000095367431640625*(double)0, t1 - t0);

			fflush(stdout);

			return EXIT_FAILURE;
		}

		printf("[%.09f] %zu bytes (%.03f MiB) were allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, sim_size, 0.00000095367431640625*(double)sim_size, t1 - t0);

		printf("[%.09f] Trying to allocate an initial %zu bytes (%.03f MiB) for the Barnes-Hut tree.\n", omp_get_wtime()-start_time, ids_size, 0.00000095367431640625*(double)ids_size);

		fflush(stdout);

		t0 = omp_get_wtime();

		ids = (int*)malloc(ids_size);

		t1 = omp_get_wtime();

		if(ids == NULL) {
			perror("Error");

			printf("[%.09f] Error: malloc() returned a NULL pointer!\n", omp_get_wtime()-start_time);

			printf("[%.09f] %zu bytes (%.03f MiB) were allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, (size_t)0, 0.00000095367431640625*(double)0, t1 - t0);

			fflush(stdout);

			return EXIT_FAILURE;
		}

		printf("[%.09f] %zu bytes (%.03f MiB) were allocated (took %.9f seconds).\n", omp_get_wtime()-start_time, ids_size, 0.00000095367431640625*(double)ids_size, t1 - t0);

		fflush(stdout);

		// partition allocated memory for the simulation
		mas = sim; // mass buffer
		pos = mas+(mas_size/sizeof(real)); // position buffer
		vel = pos+(pos_size/sizeof(real)); // velocity buffer
		acc = vel+(vel_size/sizeof(real)); // acceleration buffer
		pen = acc+(acc_size/sizeof(real)); // potential energy buffer
	}

	printf("[%.09f] Setting up initial conditions.\n", omp_get_wtime()-start_time);

	// Initialize Simulation
	{
		for(size_t i = 0; i < (N * 1); i++) mas[i] = (real)((100.0*(4.0/3.0)*PI)/(double)N);
		for(size_t i = 0; i < (N * 3); i++) pos[i] = (real)0;
		for(size_t i = 0; i < (N * 3); i++) vel[i] = (real)0;
		for(size_t i = 0; i < (N * 3); i++) acc[i] = (real)0;
		for(size_t i = 0; i < (N * 1); i++) pen[i] = (real)0;

		// randomly selected positions and velocities
		for(size_t i = 0; i < N; i++) {
			uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			//pos[3*i+0] = (real)((double)z0-0.5);
			//pos[3*i+1] = (real)((double)z1-0.5);
			//pos[3*i+2] = (real)((double)z2-0.5);

			//double ang_th = (double)TAU * (double)z0;
			//double ang_ph = acos(2.0 * (double)z1 - 1.0);
			//double ang_rh = cbrt((double)z2);
			//double sin_th = sin(ang_th);
			//double cos_th = cos(ang_th);
			//double sin_ph = sin(ang_ph);
			//double cos_ph = cos(ang_ph);
			//pos[3*i+0] = (real)(ang_rh*sin_ph*cos_th);
			//pos[3*i+1] = (real)(ang_rh*sin_ph*sin_th);
			//pos[3*i+2] = (real)(ang_rh*cos_ph);

			// normal distribution
			pos[3*i+0] = 1.000*sqrt(-2.0*log(z0))*cos(TAU*z2);
			pos[3*i+1] = 0.100*sqrt(-2.0*log(z0))*sin(TAU*z2);
			pos[3*i+2] = 1.000*sqrt(-2.0*log(z1))*cos(TAU*z3);

			//real radasdf = powf((pos[3*i+0]*pos[3*i+0])+(pos[3*i+2]*pos[3*i+2]),0.25f);
			//pos[3*i+0] /= radasdf;
			//pos[3*i+2] /= radasdf;
			//pos[3*i+0] *= 0.200f;
			//pos[3*i+1] *= 0.005f;
			//pos[3*i+2] *= 0.200f;
			//pos[2*i+0] = 0.2*sqrt(z1)*cos(TAU*z0);
			//pos[2*i+1] = 0.2*sqrt(z1)*sin(TAU*z0);
		}

		/*
		for(size_t i = 0; i < N/2; i++) {
			uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			double ang_th = (double)TAU * (double)z0;
			double ang_ph = acos(2.0 * (double)z1 - 1.0);
			double ang_rh = cbrt((double)z2);
			double sin_th = sin(ang_th);
			double cos_th = cos(ang_th);
			double sin_ph = sin(ang_ph);
			double cos_ph = cos(ang_ph);
			pos[3*i+0] = (real)(ang_rh*sin_ph*cos_th)-(real)4.0;
			pos[3*i+1] = (real)(ang_rh*sin_ph*sin_th);
			pos[3*i+2] = (real)(ang_rh*cos_ph);
		}

		for(size_t i = N/2; i < N; i++) {
			uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			double ang_th = (double)TAU * (double)z0;
			double ang_ph = acos(2.0 * (double)z1 - 1.0);
			double ang_rh = cbrt((double)z2);
			double sin_th = sin(ang_th);
			double cos_th = cos(ang_th);
			double sin_ph = sin(ang_ph);
			double cos_ph = cos(ang_ph);
			pos[3*i+0] = (real)(ang_rh*sin_ph*cos_th)+(real)4.0;
			pos[3*i+1] = (real)(ang_rh*sin_ph*sin_th);
			pos[3*i+2] = (real)(ang_rh*cos_ph);
		}*/

		for(size_t i = 0; i < N; i++) {
			vel[3*i+0] = (real)(-PI*pos[3*i+2]);
			vel[3*i+1] = (real)( 0.0*pos[3*i+1]);
			vel[3*i+2] = (real)( PI*pos[3*i+0]);
			//vel[3*i+0] = (real)0;
			//vel[3*i+1] = (real)0;
			//vel[3*i+2] = (real)0;
		}

		/*
		for(size_t i = 0; i < N/2; i++) {
			//vel[3*i+0] = (real)(-0.2*(double)pos[3*i+2]);
			//vel[3*i+1] = (real)( 0.2*(double)pos[3*i+1]);
			//vel[3*i+2] = (real)( 0.2*(double)pos[3*i+0]);
			vel[3*i+0] = (real)0.0;
			vel[3*i+1] = (real)0.0;
			//vel[3*i+2] = (real)0;
			vel[3*i+2] = (real)0;
		}
		for(size_t i = N/2; i < N; i++) {
			//vel[3*i+0] = (real)(-0.2*(double)pos[3*i+2]);
			//vel[3*i+1] = (real)( 0.2*(double)pos[3*i+1]);
			//vel[3*i+2] = (real)( 0.2*(double)pos[3*i+0]);
			vel[3*i+0] = (real)-0.0;
			vel[3*i+1] = (real)0.0;
			//vel[3*i+2] = (real)0;
			vel[3*i+2] = (real)0;
		}
		*/

		// calculate initial gravitational accelerations
		BarnesHut(tree, image, w, h, ids, mas, pos, acc, pen, N, 0);
	}

	#ifdef KEEP_CENTERED
	// center
	{
		long double avg_pos[3] = {0.0l, 0.0l, 0.0l};

		for(size_t i = 0; i < N; i++) {
			avg_pos[0] += (long double)pos[3*i+0];
			avg_pos[1] += (long double)pos[3*i+1];
			avg_pos[2] += (long double)pos[3*i+2];
		}

		avg_pos[0] /= (long double)N;
		avg_pos[1] /= (long double)N;
		avg_pos[2] /= (long double)N;

		for(size_t i = 0; i < N; i++) {
			pos[3*i+0] -= avg_pos[0];
			pos[3*i+1] -= avg_pos[1];
			pos[3*i+2] -= avg_pos[2];
		}
	}
	#endif

	printf("[%.09f] Getting ready to start the run.\n", omp_get_wtime()-start_time);

	#ifdef STARFLOOD_ENABLE_PROFILING
	FILE* diagfile = fopen("./log.csv", "w");

	fprintf(diagfile, "\"n\",\"Clear Image (ms)\",\"Render Bodies (ms)\",\"Render Quadtree (ms)\",\"Write Image (ms)\",\"First Kick (ms)\",\"Drift (ms)\",\"Update Acceleration (ms)\",\"Second Kick (ms)\"\n");

	fflush(diagfile);
	#endif

	#ifdef STARFLOOD_ENABLE_METRICS
	FILE* metricfile = fopen("./metrics.csv", "w");

	fprintf(metricfile, "\"n\",\"z\",\"kinetic energy\",\"potential energy\"\n");

	fflush(metricfile);
	#endif

	const int num_steps = FRAMES;

	float particle_lut[7*7];

	for(int i = 0; i <= 6; i++) {
		for(int j = 0; j <= 6; j++) {
			float ix = ((real)i-(real)3.)*(real)0.5;
			float iy = ((real)j-(real)3.)*(real)0.5;
			particle_lut[i*7+j] = expf(-((ix*ix)+(iy*iy)));
		}
	}

	{
		double sum = 0.;

		for(int i = 0; i < 7*7; i++) sum += particle_lut[i];

		sum = sum != 0. ? 1./sum : 1.;

		for(int i = 0; i < 7*7; i++) particle_lut[i] *= (float)sum;
	}

	printf("[%.09f] Starting run!\n", omp_get_wtime()-start_time);

	for(int step_num = 0; step_num < num_steps; step_num++) {
		printf("\r[%.09f] The run is currently on step %d/%d (z = %.09f), %.3f%% complete...", omp_get_wtime()-start_time, step_num, num_steps, TIMESTEP*(double)step_num, 100.0*((double)step_num/(double)num_steps));

		fflush(stdout);

		#ifdef STARFLOOD_ENABLE_PROFILING
		fprintf(diagfile, "%d,", step_num);
		#endif

		#ifdef STARFLOOD_ENABLE_METRICS
		{
			double ke = 0;
			double pe = 0;

			// K = (1/2)*m*v^2
			for(size_t i = 0; i < N; i++) ke += 0.5*(double)mas[i]*(double)((vel[3*i+0]*vel[3*i+0])+(vel[3*i+1]*vel[3*i+1])+(vel[3*i+2]*vel[3*i+2]));

			#ifndef METRICS_EXACT_ENERGY
			for(size_t i = 0; i < N; i++) pe += pen[i];
			#else
			for(size_t i = 0; i < N; i++) {
				for(size_t j = 0; j < N; j++) {
					if(i == j) continue;

					double disp[3] = {
						(double)pos[3*j+0] - (double)pos[3*i+0],
						(double)pos[3*j+1] - (double)pos[3*i+1],
						(double)pos[3*j+2] - (double)pos[3*i+2],
					};

					pe += G*mas[i]*mas[j]/sqrt(disp[0]*disp[0])+(disp[1]*disp[1])+(disp[2]*disp[2]); // U = G*m_i*m_j/r_ij
				}
			}
			#endif

			fprintf(metricfile, "%d,%.9f,%.9f,%.9f\n", step_num, TIMESTEP*(double)N, ke, pe);
		}
		#endif

		/*
		// reset out-of-bound particles
		for(size_t i = 0; i < N; i++) {
			if((pos[3*i+0]*pos[3*i+0])+(pos[3*i+1]*pos[3*i+1])+(pos[3*i+2]*pos[3*i+2]) > 64.0) {
				uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

				urand(&ns);

				ns += step_num;

				ns += step_num; float z0 = urand(&ns);
				
				ns += step_num; float z1 = urand(&ns);
				
				ns += step_num; float z2 = urand(&ns);
				
				ns += step_num; float z3 = urand(&ns);

				double ang_th = (double)TAU * (double)z0;
				double ang_ph = acos(2.0 * (double)z1 - 1.0);
				double ang_rh = cbrt((double)z2);
				double sin_th = sin(ang_th);
				double cos_th = cos(ang_th);
				double sin_ph = sin(ang_ph);
				double cos_ph = cos(ang_ph);
				pos[3*i+0] = (real)(8.0*(ang_rh*sin_ph*cos_th));
				pos[3*i+1] = (real)(8.0*(ang_rh*sin_ph*sin_th));
				pos[3*i+2] = (real)(8.0*(ang_rh*cos_ph));
			}
		}
		*/

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
				fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
				#endif
			}

			// render the image
			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0){
					double rms_rad = 0.;

					for(size_t i = 0; i < N; i++) {
						double rad = (pos[3*i+0]*pos[3*i+0])+(pos[3*i+1]*pos[3*i+1])+(pos[3*i+2]*pos[3*i+2]);
						rms_rad += rad;
					}

					rms_rad = 4.0*sqrt(rms_rad/(double)N);
					rms_rad = 16.0;

					for(size_t i = 0; i < N; i++) {
						// vec2 uv = (fragCoord - 0.5 * resolution) / resolution.y
						real t = 0.0166666666666666666666667*.2*(real)step_num;
						//t = 0.0;
						//real uv[2] = {pos[3 * i + 0]*cos(t)+pos[3 * i + 2]*sin(t), pos[3 * i + 1]};
						//real cam_pos[3] = {0.0f, 0.15f*(real)rms_rad, -0.9f*(real)rms_rad};
						real cam_pos[3] = {(real)0, (real)0, (real)(-1.0*(double)rms_rad)};
						real obj_dif[3] = {
						(pos[3*i+0]*cos(-t)-pos[3*i+2]*sin(-t))-cam_pos[0],
						pos[3*i+1]-cam_pos[1],
						(pos[3*i+0]*sin(-t)+pos[3*i+2]*cos(-t))-cam_pos[2]};
						/*
						real obj_dif[3] = {
						pos[3*i+0]-cam_pos[0],
						pos[3*i+1]-cam_pos[1],
						pos[3*i+2]-cam_pos[2]};
						*/
						real obj_rad = (obj_dif[0]*obj_dif[0])+(obj_dif[1]*obj_dif[1])+(obj_dif[2]*obj_dif[2]);
						if(obj_dif[2] < (real)0.0) continue;
						real uv[2] = {
						(obj_dif[0]*(real)1.)/((obj_dif[2])*((real)1.)),
						(obj_dif[1]*(real)1.)/((obj_dif[2])*((real)1.))};
						//real uv[2] = {pos[3*i+0], pos[3*i+1]};

						uv[0] *= (real)1.0; // 0.05
						uv[1] *= (real)1.0;

						int coord[2] = {
						(int)((real)h * uv[0] + (real)0.5 * (real)w),
						(int)((real)h * uv[1] + (real)0.5 * (real)h)};

						/*
						if((0 <= coord[0]) && ((size_t)coord[0] < w) && (0 <= coord[1]) && ((size_t)coord[1] < h)) {
							image[4*(w*coord[1]+coord[0])+0] += 0.05f;
							image[4*(w*coord[1]+coord[0])+1] += 0.05f;
							image[4*(w*coord[1]+coord[0])+2] += 0.05f;
							image[4*(w*coord[1]+coord[0])+3] += 0.0f;
						}
						*/
						for(int jy = -3; jy <= 3; jy++) {
							for(int jx = -3; jx <= 3; jx++) {
								int ix = jx+3;
								int iy = jy+3;
								int icoord[2] = {
								coord[0]+jx,
								coord[1]+jy};
								float br = particle_lut[iy*7+ix];

								if((0 <= icoord[0]) && ((size_t)icoord[0] < w) && (0 <= icoord[1]) && ((size_t)icoord[1] < h)) {
									image[4*(w*icoord[1]+icoord[0])+0] += 0.5*br;
									image[4*(w*icoord[1]+icoord[0])+1] += 0.5*br;
									image[4*(w*icoord[1]+icoord[0])+2] += 0.5*br;
									image[4*(w*icoord[1]+icoord[0])+3] += 0.0f;
								}
							}
						}
					}

					for(size_t i = 0; i < (w * h); i++) {
						image[4*i+3] = (real)1;
					}
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
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
				fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
				#endif
			}

			if((step_num % FRAME_INTERVAL) == 0){
				#ifdef STARFLOOD_RENDER_INTERACTIONS
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
						printf("\r[%.09f] Error: stbi_write_hdr(%s, %zu, %zu, %d, %p) returned %d (not successful)! The image may not have been written.\n", omp_get_wtime()-start_time, file_name, w, h, 4, (void*)image, 0);

						printf(
						"[%.09f] Warning: STARFLOOD_DISABLE_IMAGE_WRITE is defined,"
						"so I'm going to assume you don't want to continue wasting precious CPU time without rendering,"
						"so I'll abort the run now (at step #%d).\n", omp_get_wtime()-start_time, step_num);

						break;
					}
					#endif
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
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
			for(size_t i = 0; i < (N * 3); i++) vel[i] += (real)0.5 * dt * acc[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
			#endif
		}

		#ifdef KEEP_CENTERED
		// center
		{
			long double avg_pos[3] = {0.0l, 0.0l, 0.0l};

			for(size_t i = 0; i < N; i++) {
				avg_pos[0] += (long double)pos[3*i+0];
				avg_pos[1] += (long double)pos[3*i+1];
				avg_pos[2] += (long double)pos[3*i+2];
			}

			avg_pos[0] /= (long double)N;
			avg_pos[1] /= (long double)N;
			avg_pos[2] /= (long double)N;

			for(size_t i = 0; i < N; i++) {
				pos[3*i+0] -= avg_pos[0];
				pos[3*i+1] -= avg_pos[1];
				pos[3*i+2] -= avg_pos[2];
			}
		}
		#endif

		#ifdef HUBBLE_CONST
		for(size_t i = 0; i < (N * 3); i++) pos[i] *= (real)(1.0+0.5*(double)dt*HUBBLE_CONST);
		#endif

		// drift
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			for(size_t i = 0; i < (N * 3); i++) pos[i] += dt * vel[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
			#endif
		}

		// update acceleration
		{
			for(size_t i = 0; i < (N * 3); i++) acc[i] = (real)0;

			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			BarnesHut(tree, image, w, h, ids, mas, pos, acc, pen, N, step_num);

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.9f,", 1000.0*(t1-t0));
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
			for(size_t i = 0; i < (N * 3); i++) vel[i] += (real)0.5 * dt * acc[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.9f\n", 1000.0*(t1-t0));
			#endif
		}

		#ifdef KEEP_CENTERED
		// center
		{
			double avg_pos[3] = {0, 0, 0};

			for(size_t i = 0; i < N; i++) {
				avg_pos[0] += (double)pos[3*i+0];
				avg_pos[1] += (double)pos[3*i+1];
				avg_pos[2] += (double)pos[3*i+2];
			}

			avg_pos[0] /= (double)N;
			avg_pos[1] /= (double)N;
			avg_pos[2] /= (double)N;

			for(size_t i = 0; i < N; i++) {
				pos[3*i+0] -= avg_pos[0];
				pos[3*i+1] -= avg_pos[1];
				pos[3*i+2] -= avg_pos[2];
			}
		}
		#endif

		#ifdef HUBBLE_CONST
		for(size_t i = 0; i < (N * 3); i++) pos[i] *= (real)(1.0+0.5*(double)dt*HUBBLE_CONST);
		#endif

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

	printf("\r\x1b[2K"); // return to beginning and clear line (ANSI escape sequence used for clearing line)

	printf("[%.09f] Run complete!\n", omp_get_wtime()-start_time);

	printf("[%.09f] Cleaning up.\n", omp_get_wtime()-start_time);

	fflush(stdout);

	#ifdef STARFLOOD_ENABLE_PROFILING
	fclose(diagfile);
	#endif

	#ifdef STARFLOOD_ENABLE_METRICS
	fclose(metricfile);
	#endif

	printf("[%.09f] Attempting to free %zu bytes (%.03f MiB) used during the run.\n", omp_get_wtime()-start_time, sim_size, 0.00000095367431640625*(double)sim_size);

	fflush(stdout);

	t0 = omp_get_wtime();

	free(sim);

	t1 = omp_get_wtime();

	printf("[%.09f] %zu bytes (%.03f MiB) were freed (took %.9f seconds).\n", omp_get_wtime()-start_time, sim_size, 0.00000095367431640625*(double)sim_size, t1 - t0);

	printf("[%.09f] Attempting to free %zu bytes (%.03f MiB) used during rendering.\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size);

	fflush(stdout);

	t0 = omp_get_wtime();

	free(image);

	t1 = omp_get_wtime();

	printf("[%.09f] %zu bytes (%.03f MiB) were freed (took %.9f seconds).\n", omp_get_wtime()-start_time, image_size, 0.00000095367431640625*(double)image_size, t1 - t0);

	printf("[%.09f] Execution complete. Thank you for using Starflood!\n", omp_get_wtime()-start_time);

	fflush(stdout);

	return EXIT_SUCCESS;
}
