#include <simulation/simulation.hpp>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <utilities/graphics.hpp>
#include <utilities/common.h>

#include <cstdio>
#include <cstdlib>

void sim::run() {

	const real dt = TIMESTEP; // Simulation Timestep


	auto mas = sim_particles.mas;
	auto pos = sim_particles.pos;
	auto vel = sim_particles.vel;
	auto acc = sim_particles.acc;
	auto pen = sim_particles.pen;
	auto ids = sim_particles.ids;

	printf("[%.09f] Getting ready to start the run.\n", omp_get_wtime()-start_time);

	#ifdef STARFLOOD_ENABLE_PROFILING
	diagfile = fopen("./log.csv", "w");

	fprintf(diagfile, "\"n\",\"Clear Image (ms)\",\"Render Bodies (ms)\",\"Render Quadtree (ms)\",\"Write Image (ms)\",\"First Kick (ms)\",\"Drift (ms)\",\"Update Acceleration (ms)\",\"Second Kick (ms)\"\n");

	fflush(diagfile);
	#endif

	#ifdef STARFLOOD_ENABLE_STATISTICS
	statistics_file = fopen("./statistics.csv", "w");

	fprintf(statistics_file, "\"n\"\t\"z\"\t\"kinetic energy\"\t\"potential energy\"\n");

	fflush(statistics_file);
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
		printf("\r[%.09f] The run is currently on step %d/%d (z = %.03f), %.3f%% complete...", omp_get_wtime()-start_time, step_num, num_steps, (TIMESTEP*(double)num_steps)-(TIMESTEP*(double)step_num), 100.0*((double)step_num/(double)num_steps));

		fflush(stdout);

		#ifdef STARFLOOD_ENABLE_PROFILING
		fprintf(diagfile, "%d,", step_num);
		#endif

		#ifdef STARFLOOD_ENABLE_STATISTICS
		{
			double ke = 0;
			double pe = 0;

			// K = (1/2)*m*v^2
			for(size_t i = 0; i < N; i++) ke += 0.5*(double)mas[i]*(double)((vel[3*i+0]*vel[3*i+0])+(vel[3*i+1]*vel[3*i+1])+(vel[3*i+2]*vel[3*i+2]));

			#ifndef STATISTICS_EXACT_ENERGY
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

					pe += G*mas[i]*mas[j]/sqrt((disp[0]*disp[0])+(disp[1]*disp[1])+(disp[2]*disp[2])+SOFTENING_PARAMETER); // U = G*m_i*m_j/r_ij
				}
			}
			#endif

			fprintf(statistics_file, "%d\t%.9f\t%.9f\t%.9f\n", step_num, (TIMESTEP*(double)num_steps)-(TIMESTEP*(double)step_num), ke, pe);
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

			double rms_rad = 0.; // root mean square radius

			for(size_t i = 0; i < N; i++) rms_rad += (pos[3*i+0]*pos[3*i+0])+(pos[3*i+1]*pos[3*i+1])+(pos[3*i+2]*pos[3*i+2]);

			rms_rad = 4.0*sqrt(rms_rad/(double)N);
			//rms_rad = 2.0;//16.0;

			// camera position
			//real cam_pos[3] = {0.0f, 0.15f*(real)rms_rad, -0.9f*(real)rms_rad};
			real cam_pos[3] = {(real)0, (real)(0.25*(double)rms_rad), (real)(-1.0*(double)rms_rad)};

			real cam_dir[3] = {(real)0.0, (real)0.0, (real)1.0};

			// camera time
			real time = 0.0166666666666666666666667*0.*(real)step_num;

			real focal_length = 1;

			// render the image
			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0) {
					for(size_t i = 0; i < N; i++) {

						//real obj_dif[3] = {
						//(vert[3*i+0]*cos(-t)-vert[3*i+2]*sin(-t))-cam_pos[0],
						//vert[3*i+1]-cam_pos[1],
						//(vert[3*i+0]*sin(-t)+vert[3*i+2]*cos(-t))-cam_pos[2]};

						//real t = 0.0166666666666666666666667*0.*(real)step_num;
						//t = 0.0;
						//real uv[2] = {pos[3 * i + 0]*cos(t)+pos[3 * i + 2]*sin(t), pos[3 * i + 1]};

						real vout[3];
						real vert[3] = {pos[3*i+0], pos[3*i+1], pos[3*i+2]};

						projection(vout, vert, cam_pos, cam_dir, focal_length, time);

						if(vout[2] < (real)0.0) continue; // particle is behind camera

						// inverse of vec2 uv = (fragCoord - 0.5 * resolution) / resolution.y
						int coord[2] = {
							(int)((real)h * vout[0] + (real)0.5 * (real)w),
							(int)((real)h * vout[1] + (real)0.5 * (real)h)
						};

						/*
						// single pixel
						if((0 <= coord[0]) && ((size_t)coord[0] < w) && (0 <= coord[1]) && ((size_t)coord[1] < h)) {
							image[4*(w*coord[1]+coord[0])+0] += 0.05f;
							image[4*(w*coord[1]+coord[0])+1] += 0.05f;
							image[4*(w*coord[1]+coord[0])+2] += 0.05f;
							image[4*(w*coord[1]+coord[0])+3] += 0.0f;
						}
						*/
						// gaussian
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

					// clear alpha (opaque)
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
					{
					{
						real p0[3] = {(real)(-1.0), (real)(-1.0), (real)(-1.0)};
						real p1[3] = {(real)(-1.0), (real)(-1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)( 1.0), (real)(-1.0), (real)(-1.0)};
						real p1[3] = {(real)( 1.0), (real)(-1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)(-1.0), (real)( 1.0), (real)(-1.0)};
						real p1[3] = {(real)(-1.0), (real)( 1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)( 1.0), (real)( 1.0), (real)(-1.0)};
						real p1[3] = {(real)( 1.0), (real)( 1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}

					{
						real p0[3] = {(real)(-1.0), (real)(-1.0), (real)(-1.0)};
						real p1[3] = {(real)( 1.0), (real)(-1.0), (real)(-1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)(-1.0), (real)(-1.0), (real)(-1.0)};
						real p1[3] = {(real)(-1.0), (real)( 1.0), (real)(-1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)( 1.0), (real)( 1.0), (real)(-1.0)};
						real p1[3] = {(real)( 1.0), (real)(-1.0), (real)(-1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)( 1.0), (real)( 1.0), (real)(-1.0)};
						real p1[3] = {(real)(-1.0), (real)( 1.0), (real)(-1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}

					{
						real p0[3] = {(real)(-1.0), (real)(-1.0), (real)( 1.0)};
						real p1[3] = {(real)( 1.0), (real)(-1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)(-1.0), (real)(-1.0), (real)( 1.0)};
						real p1[3] = {(real)(-1.0), (real)( 1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)( 1.0), (real)( 1.0), (real)( 1.0)};
						real p1[3] = {(real)( 1.0), (real)(-1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					{
						real p0[3] = {(real)( 1.0), (real)( 1.0), (real)( 1.0)};
						real p1[3] = {(real)(-1.0), (real)( 1.0), (real)( 1.0)};
						drawLine3D(image, w, h, 0.250f, 1.000f, 0.250f, 1.000f, cam_pos, cam_dir, focal_length, time, p0, p1);
					}
					}

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
}
