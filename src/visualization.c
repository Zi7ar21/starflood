// Needed for pthread.h and posix_memalign()
#define _POSIX_C_SOURCE 200809L

#include "visualization.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "common.h"
#include "config.h"
#include "log.h"
#include "rng.h"
#include "simulation.h"
#include "solver.h"
#include "types.h"
#include "timing.h"

#ifdef ENABLE_VIS_IO_THREAD
#include <pthread.h>
extern pthread_t vis_io_thread;
extern pthread_mutex_t vis_io_mutex;

void* dummy_function(void* arg) {
	return arg;
}
#endif

#ifdef LOG_TIMINGS_VIS_DRAW
log_t log_timings_vis_draw;
#endif

int visualization_init(vis_t* restrict visualization, unsigned int w, unsigned int h) {
	#ifdef ENABLE_VIS_IO_THREAD
	if( 0 != pthread_create(&vis_io_thread, NULL, dummy_function, NULL) ) {
		fprintf(stderr, "%s error: pthread_create(&vis_io_thread, NULL, dummy_function, NULL) ", "visualization_init()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	if( 0 != pthread_mutex_init(&vis_io_mutex, NULL) ) {
		fprintf(stderr, "%s error: pthread_mutex_init(&vis_io_mutex, NULL) ", "visualization_init()");
		perror("failed");
		pthread_join(vis_io_thread, NULL);
		return STARFLOOD_FAILURE;
	}
	#endif

	TIMING_INIT();

	#ifdef LOG_TIMINGS_VIS_DRAW
	if( STARFLOOD_SUCCESS != log_init(&log_timings_vis_draw, OUTPUT_DIR "/" LOG_TIMINGS_VIS_DRAW) ) {
		return STARFLOOD_FAILURE;
	}

	fprintf(log_timings_vis_draw.file, "%s,%s,%s,%s,%s\n", "step_number", "clear_atomic_buffer", "calc_matMVP", "rasterize_atomic", "finalize");

	log_sync(&log_timings_vis_draw);
	#endif

	vis_t vis = *visualization;

	void* mem = NULL;

	i32* atomic_buffer = (i32*)NULL;
	f32* render_buffer = (f32*)NULL;

	#if (0 >= VIS_FILE_FORMAT)
	f32* binary_buffer = (f32*)NULL;
	#else
	unsigned char* binary_buffer = (unsigned char*)NULL;
	#endif

	printf("Visualization Memory Addresses:\n");

	size_t atomic_buffer_length = (size_t)4u * (size_t)w * (size_t)h; //  32-bit        integer, RGBA, for atomic accumulation
	size_t render_buffer_length = (size_t)4u * (size_t)w * (size_t)h; //  32-bit floating-point, RGBA, for finishing the render and image post-processing (alpha channel is currently unused, only for alignment)
	size_t binary_buffer_length = (size_t)3u * (size_t)w * (size_t)h; // image format-dependent, RGB , for serializing the finalized image

	#if (0 >= VIS_FILE_FORMAT)
	size_t mem_size = (sizeof(i32)*atomic_buffer_length)+(sizeof(i32)*render_buffer_length)+(sizeof(f32)*binary_buffer_length);
	#else
	size_t mem_size = (sizeof(i32)*atomic_buffer_length)+(sizeof(i32)*render_buffer_length)+(sizeof(unsigned char)*binary_buffer_length);
	#endif

	TIMING_START();

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
	#endif

	TIMING_STOP();

	#ifdef STARFLOOD_ALIGNMENT
	TIMING_PRINT("visualization_init()", "posix_memalign()");
	#else
	TIMING_PRINT("visualization_init()", "malloc()");
	#endif

	if(NULL == mem) {
		#ifdef STARFLOOD_ALIGNMENT
		fprintf(stderr, "%s error: mem is NULL after posix_memalign(&mem, %zu, %zu", "visualization_init()", (size_t)STARFLOOD_ALIGNMENT, mem_size);
		#else
		fprintf(stderr, "%s error: mem is NULL after malloc(%zu", "visualization_init()", mem_size);
		#endif

		perror(")");

		return EXIT_FAILURE;
	}

	printf("  mem: %p (%zu bytes)\n", mem, mem_size);

	TIMING_START();

	memset(mem, 0, mem_size);

	TIMING_STOP();
	TIMING_PRINT("visualization_init()", "memset()");

	atomic_buffer = (i32*)mem;
	render_buffer = (f32*)&atomic_buffer[atomic_buffer_length];
	#if (0 >= VIS_FILE_FORMAT)
	binary_buffer = (f32*)&render_buffer[render_buffer_length];
	#else
	binary_buffer = (unsigned char*)&render_buffer[render_buffer_length];
	#endif

	printf("  atomic_buffer: %p (%zu bytes)\n", (void*)atomic_buffer, sizeof(i32)*atomic_buffer_length);
	printf("  render_buffer: %p (%zu bytes)\n", (void*)render_buffer, sizeof(f32)*render_buffer_length);
	#if (0 >= VIS_FILE_FORMAT)
	printf("  binary_buffer: %p (%zu bytes)\n", (void*)binary_buffer, sizeof(f32)*binary_buffer_length);
	#else
	printf("  binary_buffer: %p (%zu bytes)\n", (void*)binary_buffer, sizeof(unsigned char)*binary_buffer_length);
	#endif
	printf("\n");

	TIMING_START();

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = (i32)0;
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_init()", "initialize atomic_buffer");
	TIMING_START();

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		render_buffer[i] = (f32)0.000;
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_init()", "initialize render_buffer");

	TIMING_START();
	for(unsigned int i = 0u; i < 3u * w * h; i++) {
		#if (0 >= VIS_FILE_FORMAT)
		binary_buffer[i] = (f32)0.000;
		#else
		binary_buffer[i] = (unsigned char)0u;
		#endif
	}
	TIMING_STOP();
	TIMING_PRINT("visualization_init()", "initialize binary_buffer");

	vis.w = w;
	vis.h = h;
	vis.mem = mem;
	vis.atomic_buffer = atomic_buffer;
	vis.render_buffer = render_buffer;
	vis.binary_buffer = binary_buffer;

	*visualization = vis;

	return STARFLOOD_SUCCESS;
}

int visualization_free(vis_t* restrict visualization) {
	#ifdef ENABLE_VIS_IO_THREAD
	pthread_join(vis_io_thread, NULL);
	pthread_mutex_destroy(&vis_io_mutex);
	#endif

	TIMING_INIT();

	#ifdef LOG_TIMINGS_VIS_DRAW
	if( STARFLOOD_SUCCESS != log_free(&log_timings_vis_draw) ) {
		fprintf(stderr, "%s error: %s failed!\n", "visualization_free()", "log_free(&log_timings_vis_draw)");
	}
	#endif

	vis_t vis = *visualization;

	TIMING_START();

	free(vis.mem);

	TIMING_STOP();
	TIMING_PRINT("visualization_free()", "free()");

	vis.mem = NULL;

	vis.atomic_buffer = (i32*)NULL;
	vis.render_buffer = (f32*)NULL;

	#if (0 >= VIS_FILE_FORMAT)
	vis.binary_buffer = (f32*)NULL;
	#else
	vis.binary_buffer = (unsigned char*)NULL;
	#endif

	*visualization = vis;

	return STARFLOOD_SUCCESS;
}

int visualization_draw(const vis_t* restrict visualization, const sim_t* restrict sim_ptr) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	vis_t vis = *visualization;

	unsigned int w = vis.w;
	unsigned int h = vis.h;

	i32* atomic_buffer = vis.atomic_buffer;
	f32* render_buffer = vis.render_buffer;

	unsigned int step_number = sim.step_number;
	unsigned int           N = sim.N;

	real* pos = sim_find(&sim, SIM_POS);
	real* vel = sim_find(&sim, SIM_VEL);
	real* mas = sim_find(&sim, SIM_MAS);
	real* pot = sim_find(&sim, SIM_POT);

	double average_radius = 1.0;

	{
		average_radius = 0.0;

		double sum_com = 0.0;

		for(unsigned int i = 0u; i < N; i++) {
			double pos_i[3] = {
				(double)pos[3u*i+0u],
				(double)pos[3u*i+1u],
				(double)pos[3u*i+2u]
			};

			double r2 = (pos_i[0]*pos_i[0])+(pos_i[1]*pos_i[1])+(pos_i[2]*pos_i[2]);
			double r1 = sqrt(r2);

			{
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				double y = r1 - sum_com;
				volatile double t = average_radius + y;
				volatile double z = t - average_radius;
				sum_com = z - y;
				average_radius = t;

				// Naïve summation
				//average_radius += r1;
			}
		}

		average_radius *= 1.0 / (double)N;

		printf("average_radius: % .15f\n", average_radius);
	}

	double pot_min = (double)pot[0u];
	double pot_max = (double)pot[0u];

	for(unsigned int i = 1u; i < N; i++) {
		pot_min = fmin(pot_min, (double)pot[i]);
		pot_max = fmax(pot_max, (double)pot[i]);
	}

	#ifdef LOG_TIMINGS_VIS_DRAW
	fprintf(log_timings_vis_draw.file, "%u", step_number);
	#endif

	#ifdef OUTPUT_INTERVAL
	real time = (real)( (double)TIMESTEP_SIZE * (double)step_number );
	#else
	real time = (real)( (double)TIMESTEP_SIZE * (double)step_number );
	#endif

	TIMING_START();

	#ifdef VISUALIZATION_RASTERIZATION
	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = (i32)0;
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_draw()", "clear_atomic_buffer");
	#ifdef LOG_TIMINGS_VIS_DRAW
	LOG_TIMING(log_timings_vis_draw);
	#endif
	TIMING_START();

	// model matrix
	real matM[16];

	// model matrix = rotation matrix
	{
		// I feel genuinely sorry for the compiler
		// that has to optimize all of this lol

		real alp = (real)0.0; // roll
		real bet = (real)0.0; // pitch
		real gam = (real)0.0; // yaw

		//bet = (real)TAU * (real)0.0000;
		//bet = (real)TAU * (real)0.0625;
		//bet = (real)TAU * (real)0.1250;
		//bet = (real)TAU * (real)0.2500;
		//bet = (real)TAU * (real)1.000e-3 * time;

		//gam = (real)TAU * (real)0.0000;
		//gam = (real)TAU * (real)0.0025;
		//gam = (real)TAU * (real)0.0625;
		//gam = (real)TAU * (real)0.1250;
		gam = (real)TAU * (real)0.2500;
		//gam = (real)TAU * (real)0.0625 * cos( (real)1.000e-2 * time * (real)TAU );
		//gam = (real)TAU * (real)1.000e-3 * time;

		real cos_alp = real_cos(alp), sin_alp = real_sin(alp);
		real cos_bet = real_cos(bet), sin_bet = real_sin(bet);
		real cos_gam = real_cos(gam), sin_gam = real_sin(gam);

		real rot_mat[16] = {
			(real)(cos_alp*cos_bet), (real)(cos_alp*sin_bet*sin_gam-sin_alp*cos_gam), (real)(cos_alp*sin_bet*cos_gam+sin_alp*sin_gam), (real)(0.0),
			(real)(sin_alp*cos_bet), (real)(sin_alp*sin_bet*sin_gam+cos_alp*cos_gam), (real)(sin_alp*sin_bet*cos_gam-cos_alp*sin_gam), (real)(0.0),
			(real)(       -sin_bet), (real)(                        cos_bet*sin_gam), (real)(                        cos_bet*cos_gam), (real)(0.0),
			(real)(            0.0), (real)(                                    0.0), (real)(                                    0.0), (real)(1.0)
		};

		for(int i = 0; i < 16; i++) {
			matM[i] = rot_mat[i];
		}
	}

	// view matrix
	real matV[16] = {
		(real)( 1.000), (real)( 0.000), (real)( 0.000), (real)( 0.000),
		(real)( 0.000), (real)( 1.000), (real)( 0.000), (real)( 0.000),
		(real)( 0.000), (real)( 0.000), (real)( 1.000), (real)( 0.000),
		(real)( 0.000), (real)( 0.000), (real)( 0.000), (real)( 1.000)
	};

	#ifdef ORTHO_SCALE
	real ortho_scale = (real)exp2(ORTHO_SCALE);
	#else
	real ortho_scale = (real)1.0;
	#endif

	#ifdef AUTOSCALING
	ortho_scale *= (real)(1.0 / average_radius);
	#endif

	// projection matrix
	real matP[16] = {
		(real)(ortho_scale), (real)(        0.0), (real)(        0.0), (real)(        0.0),
		(real)(        0.0), (real)(ortho_scale), (real)(        0.0), (real)(        0.0),
		(real)(        0.0), (real)(        0.0), (real)(ortho_scale), (real)(        0.0),
		(real)(        0.0), (real)(        0.0), (real)(        0.0), (real)(        1.0)
	};

	real matMVP[16];

	// C = A * B
	{
		real A[16], B[16], C[16];

		for(int i = 0; i < 16; i++) {
			A[i] = matV[i];
		}

		for(int i = 0; i < 16; i++) {
			B[i] = matM[i];
		}

		// matrix multiplication
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				C[4*i+j] =
				(A[4*i+3]*B[4*3+j])+
				(A[4*i+2]*B[4*2+j])+
				(A[4*i+1]*B[4*1+j])+
				(A[4*i+0]*B[4*0+j]);
			}
		}

		for(int i = 0; i < 16; i++) {
			matMVP[i] = C[i];
		}
	}

	// C = A * B
	{
		real A[16], B[16], C[16];

		for(int i = 0; i < 16; i++) {
			A[i] = matP[i];
		}

		for(int i = 0; i < 16; i++) {
			B[i] = matMVP[i];
		}

		// matrix multiplication
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				C[4*i+j] =
				(A[4*i+3]*B[4*3+j])+
				(A[4*i+2]*B[4*2+j])+
				(A[4*i+1]*B[4*1+j])+
				(A[4*i+0]*B[4*0+j]);
			}
		}

		for(int i = 0; i < 16; i++) {
			matMVP[i] = C[i];
		}
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_draw()", "calc_matMVP");
	#ifdef LOG_TIMINGS_VIS_DRAW
	LOG_TIMING(log_timings_vis_draw);
	#endif
	TIMING_START();

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_VIS
		#pragma omp target teams distribute parallel for map(tofrom: atomic_buffer[:4u*w*h]) map(to: pos[:3u*N], vel[:3u*N], pot[:N])
		#else
		#pragma omp parallel for schedule(dynamic, 256)
		#endif
	#endif
	for(unsigned int idx = 0u; idx < N; idx++) {
		// https://cs418.cs.illinois.edu/website/text/math2.html

		// the fourth component is for homogenous coordinates
		// https://en.wikipedia.org/wiki/Homogeneous_coordinates
		real p[4] = {
			(real)(pos[3u*idx+0u]),
			(real)(pos[3u*idx+1u]),
			(real)(pos[3u*idx+2u]),
			(real)(           1.0)
		};

		real v[3] = {
			vel[3u*idx+0u],
			vel[3u*idx+1u],
			vel[3u*idx+2u]
		};

		real potential = pot[idx];

		// Apply (Model, View, Projection) Matrix
		// y = A * x
		{
			real A[16], x[4], y[4];

			for(int i = 0; i < 16; i++) {
				A[i] = matMVP[i];
			}

			for(int i = 0; i < 4; i++) {
				x[i] = p[i];
			}

			for(int i = 0; i < 4; i++) {
				y[i] = (A[4*i+0]*x[0])+(A[4*i+1]*x[1])+(A[4*i+2]*x[2])+(A[4*i+3]*x[3]);
			}

			for(int i = 0; i < 4; i++) {
				p[i] = y[i];
			}
		}

		real palette_phase = (real)1.000e-1 * potential;

		real color[3] = {
			(real)0.5*real_cos( (real)TAU * palette_phase - (real)( TAU * (0.0/3.0) ) )+(real)0.5,
			(real)0.5*real_cos( (real)TAU * palette_phase - (real)( TAU * (1.0/3.0) ) )+(real)0.5,
			(real)0.5*real_cos( (real)TAU * palette_phase - (real)( TAU * (2.0/3.0) ) )+(real)0.5
		};

		for(unsigned int samp_number = 0u; samp_number < (unsigned int)SPATIAL_SAMPLES; samp_number++) {
			i32 accumulation_color[4] = {
				(i32)( (real)100.0 * color[0] ),
				(i32)( (real)100.0 * color[1] ),
				(i32)( (real)100.0 * color[2] ),
				(i32)1
			};

			real p_i[3] = {
				p[0],
				p[1],
				p[2]
			};

			real time_offset = (real)0.0;

			real sample_offset[2] = {(real)0.0, (real)0.0};

			u32 s[4] = {
				(u32)0xB79ABC95u + (u32)samp_number,
				(u32)0xE0DA3F84u + (u32)idx,
				(u32)0xAB75F07Bu + (u32)step_number,
				(u32)0xCF52CA01u
			};

			pcg4d(s);
			//pcg4d(s); // second round for better statistical quality

			{
				real r[4] = {
					(real)INV_PCG32_MAX * (real)s[0],
					(real)INV_PCG32_MAX * (real)s[1],
					(real)INV_PCG32_MAX * (real)s[2],
					(real)INV_PCG32_MAX * (real)s[3]
				};

				// Box-Muller Transform
				// https://en.wikipedia.org/wiki/Box–Muller_transform
				real n[2] = {
					real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_cos( (real)TAU * r[3] ),
					real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_sin( (real)TAU * r[3] )
				};

				// Spatial anti-aliasing (gaussian window)
				{
					const real filter_sigma = (real)0.750;

					sample_offset[0] = filter_sigma * n[0];
					sample_offset[1] = filter_sigma * n[1];
				}

				{
					#ifdef SHUTTER_SPEED
					time_offset += (real)SHUTTER_SPEED * (r[0] - (real)0.5); // motion blur
					#endif

					#ifdef OUTPUT_INTERVAL
					time_offset *= (real)OUTPUT_INTERVAL * (real)TIMESTEP_SIZE;
					#else
					time_offset *= (real)TIMESTEP_SIZE;
					#endif

					// predicted position using velocity vector
					for(int i = 0; i < 3; i++) {
						p_i[i] += time_offset * v[i];
					}
				}
			}

			/*
			// rotate2d
			{
				real x[2], y[2];

				x[0] = p_i[1]; // Y
				x[1] = p_i[2]; // Z

				//const real theta = (real)TAU * (real)0.0000; // side-on
				const real theta = (real)TAU * (real)0.0625; // slightly tilted
				//const real theta = (real)TAU * (real)0.1250; // edge-on
				//const real theta = (real)TAU * (real)0.2500; // top-down

				real cos_theta = real_cos(theta), sin_theta = real_sin(theta);

				y[0] = (cos_theta*x[0])-(sin_theta*x[1]);
				y[1] = (sin_theta*x[0])+(cos_theta*x[1]);

				p_i[1] = y[0];
				p_i[2] = y[1];
			}

			// rotate2d
			{
				real x[2], y[2];

				x[0] = p_i[0]; // X
				x[1] = p_i[2]; // Z

				//const real theta = (real)TAU * (real)0.0000;
				//const real theta = (real)TAU * (real)0.0625;
				//const real theta = (real)TAU * (real)0.1250;
				//const real theta = (real)TAU * (real)0.2500;
				real theta = (real)TAU * 0.001 * (time + time_offset);

				real cos_theta = real_cos(theta), sin_theta = real_sin(theta);

				y[0] = (cos_theta*x[0])-(sin_theta*x[1]);
				y[1] = (sin_theta*x[0])+(cos_theta*x[1]);

				p_i[0] = y[0];
				p_i[2] = y[1];
			}
			*/

			int coord[2] = {
				(int)(sample_offset[0] - (real)0.5 + (real)0.5 * (real)h * p_i[0] + (real)0.5 * (real)w),
				(int)(sample_offset[1] - (real)0.5 + (real)0.5 * (real)h * p_i[1] + (real)0.5 * (real)h)
			};

			if( 0 <= coord[0] && coord[0] < (int)w && 0 <= coord[1] && coord[1] < (int)h ) {
				unsigned int pixel_index = w * (unsigned int)coord[1] + (unsigned int)coord[0];

				for(unsigned int i = 0u; i < 4u; i++) {
					#ifdef _OPENMP
					#pragma omp atomic update
					#endif
					atomic_buffer[4u*pixel_index+i] += accumulation_color[i];
				}
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_draw()", "rasterize_atomic");
	#ifdef LOG_TIMINGS_VIS_DRAW
	LOG_TIMING(log_timings_vis_draw);
	#endif
	TIMING_START();

	// floating-point scaling factor for integer atomic_buffer
	double pixel_value_scale = 0.01 / (double)SPATIAL_SAMPLES;

	#ifdef EXPOSURE
	pixel_value_scale *= exp2(EXPOSURE);
	#endif

	TIMING_START();

		// read the atomic buffer and finish rendering the visualization
	for(unsigned int i = 0u; i < w * h; i++) {
		double pixel[4] = {
			(double)0.250,
			(double)0.250,
			(double)0.250,
			(double)1.000
		};

		i32 accum[4];

		for(unsigned int j = 0u; j < 4u; j++) {
			i32 val = (i32)0;

			#ifdef _OPENMP
			#pragma omp atomic read
			#endif
			val = atomic_buffer[4u*i+j];

			accum[j] = val;
		}

		for(int j = 0; j < 4; j++) {
			pixel[j] = pixel_value_scale * (double)accum[j];
		}

		#if (0 < VIS_FILE_FORMAT)
		// tonemapping
		for(int j = 0; j < 4; j++) {
			pixel[j] = tanh(pixel[j]);
		}
		#endif

		// set alpha channel to 1.000
		pixel[3] = (double)1.000;

		for(unsigned int j = 0u; j < 4u; j++) {
			render_buffer[4u*i+j] = (f32)pixel[j];
		}
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_draw()", "finalize");
	#endif

	#if 0
	real uv_normalization_factor = (real)(0.5 / (double)h);

	#ifdef ORTHO_SCALE
	const real ortho_scale = (real)exp2(-ORTHO_SCALE);
	#else
	const real ortho_scale = (real)1.0;
	#endif

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_VIS
		#pragma omp target teams distribute parallel for collapse(2) map(tofrom: render_buffer[:4u*w*h]) map(to: pos[:3u*N], mas[:N])
		#else
		#pragma omp parallel for collapse(2) schedule(dynamic, 256)
		#endif
		//#pragma omp tile sizes(16,16)
	#endif
	for(unsigned int idx = 0u; idx < w; idx++) {
		for(unsigned int idy = 0u; idy < h; idy++) {
			real pixel_color[3] = {(real)0.0, (real)0.0, (real)0.0};

			for(unsigned int samp_number = 0u; samp_number < (unsigned int)SPATIAL_SAMPLES; samp_number++) {
				u32 s[4] = {
					(u32)0xB79ABC95u + (u32)samp_number,
					(u32)0xE0DA3F84u + (u32)(w*idy+idx),
					(u32)0xAB75F07Bu + (u32)step_number,
					(u32)0xCF52CA01u
				};

				pcg4d(s);
				//pcg4d(s); // second round for better statistical quality

				real r[4] = {
					(real)INV_PCG32_MAX * (real)s[0],
					(real)INV_PCG32_MAX * (real)s[1],
					(real)INV_PCG32_MAX * (real)s[2],
					(real)INV_PCG32_MAX * (real)s[3]
				};

				// Box-Muller Transform
				// https://en.wikipedia.org/wiki/Box–Muller_transform
				real n[2] = {
					real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_cos( (real)TAU * r[3] ),
					real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_sin( (real)TAU * r[3] )
				};

				// Spatial anti-aliasing (gaussian window)
				const real filter_sigma = (real)0.750;

				real uv[2] = {
					uv_normalization_factor * (filter_sigma * n[0] + (real)0.5 + (real)idx - (real)0.5 * (real)w ),
					uv_normalization_factor * (filter_sigma * n[1] + (real)0.5 + (real)idy - (real)0.5 * (real)h )
				};

				uv[0] *= ortho_scale;
				uv[1] *= ortho_scale;

				/*
				real probe_pos[3] = {
					(real)( uv[0]),
					(real)( uv[1]),
					(real)( 0.000)
				};
				*/

				real probe_pos[3] = {
					(real)( uv[0]),
					(real)( 0.000),
					(real)( uv[1])
				};

				/*
				real probe_pos[3] = {
					(real)(real_cos( (real)TAU * r[1]) * uv[0]),
					(real)(                              uv[1]),
					(real)(real_sin( (real)TAU * r[1]) * uv[0])
				};
				*/

				//real potential = probe_potential(mas, pos, probe_pos, N, step_number);
				real potential = (real)0.0;

				{
	#ifdef PROBE_DECIMATION
	unsigned int i_length = N / (unsigned int)PROBE_DECIMATION;
	unsigned int i_offset = (step_number % (unsigned int)PROBE_DECIMATION) * i_length;
	#endif

	real acc_sum[3] = {(real)0.0, (real)0.0, (real)0.0};
	real pot_sum = (real)0.0;
	#ifdef ENABLE_KAHAN_SUMMATION
	real acc_com[3] = {(real)0.0, (real)0.0, (real)0.0};
	real pot_com = (real)0.0;
	#endif

	#ifdef PROBE_DECIMATION
	for(unsigned int i = i_offset; i < (i_offset+i_length); i++) {
	#else
	for(unsigned int i = 0u; i < N; i++) {
	#endif
		real m_i = mas[i];

		real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		real r_ij[3] = {
			probe_pos[0] - r_i[0],
			probe_pos[1] - r_i[1],
			probe_pos[2] - r_i[2]
		};

		real r2 = (r_ij[0u]*r_ij[0u])+(r_ij[1u]*r_ij[1u])+(r_ij[2u]*r_ij[2u]);

		#ifdef EPSILON
		real inv_r2 = (real)1.0 /          ( r2 + (real)(EPSILON*EPSILON) );
		real inv_r1 = (real)1.0 / real_sqrt( r2 + (real)(EPSILON*EPSILON) );
		#else
		real inv_r2 = (real)0.0 < r2 ? (real)1.0 /          (r2) : (real)0.0;
		real inv_r1 = (real)0.0 < r2 ? (real)1.0 / real_sqrt(r2) : (real)0.0;
		#endif

		// gravitational potential of body j
		// (G and m_i are taken into account later)
		real pot_i = m_i * inv_r1;

		real F[3] = {
			pot_i * r_ij[0u] * inv_r2,
			pot_i * r_ij[1u] * inv_r2,
			pot_i * r_ij[2u] * inv_r2
		};

		for(unsigned int k = 0u; k < 3u; k++) {
			#ifdef ENABLE_KAHAN_SUMMATION
			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			real y = F[k] - acc_com[k];
			volatile real t = acc_sum[k] + y;
			volatile real z = t - acc_sum[k];
			acc_com[k] = z - y;
			acc_sum[k] = t;
			#else
			// Naïve summation
			acc_sum[k] += F[k];
			#endif
		}

		{
			#ifdef ENABLE_KAHAN_SUMMATION
			// Kahan summation
			// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
			real y = pot_i - pot_com;
			volatile real t = pot_sum + y;
			volatile real z = t - pot_sum;
			pot_com = z - y;
			pot_sum = t;
			#else
			// Naïve summation
			pot_sum += pot_j;
			#endif
		}
	}

	#ifdef PROBE_DECIMATION
	//acc[3u*i+0u] = (real)PROBE_DECIMATION * (real)(-G) * m_i * acc_sum[0u];
	//acc[3u*i+1u] = (real)PROBE_DECIMATION * (real)(-G) * m_i * acc_sum[1u];
	//acc[3u*i+2u] = (real)PROBE_DECIMATION * (real)(-G) * m_i * acc_sum[2u];

	//return (real)PROBE_DECIMATION * (real)(-G) * pot_sum;
	potential = (real)PROBE_DECIMATION * (real)(-G) * pot_sum;
	#else
	//acc[3u*i+0u] = (real)(-G) * m_i * acc_sum[0u];
	//acc[3u*i+1u] = (real)(-G) * m_i * acc_sum[1u];
	//acc[3u*i+2u] = (real)(-G) * m_i * acc_sum[2u];

	//return (real)(-G) * pot_sum;
	potential = (real)(-G) * pot_sum;
	#endif
				}

				real color_phase = (real)1.000e1 * potential;

				real color[3] = {
					(real)0.5*real_cos( (real)TAU * color_phase - (real)( TAU * (0.0/3.0) ) )+(real)0.5,
					(real)0.5*real_cos( (real)TAU * color_phase - (real)( TAU * (1.0/3.0) ) )+(real)0.5,
					(real)0.5*real_cos( (real)TAU * color_phase - (real)( TAU * (2.0/3.0) ) )+(real)0.5
				};

				for(int i = 0; i < 3; i++) {
					//color[i] *= log10( -potential + ( 1.0 / pow(10.0, 1.0) ) ) + 1.0;
					//color[i] *= -1.000e-1 * potential;
					color[i] = -1.000e-1 * potential;
				}

				for(int i = 0; i < 3; i++) {
					pixel_color[i] += color[i];
					//pixel_color[i] = fmaxf(pixel_color[i], color[i]);
				}
			}

			real scale_factor = (real)(1.0 / (double)SPATIAL_SAMPLES);

			for(int i = 0; i < 3; i++) {
				pixel_color[i] *= scale_factor;
			}

			render_buffer[4u*(w*idy+idx)+0u] = (f32)(pixel_color[0u]);
			render_buffer[4u*(w*idy+idx)+1u] = (f32)(pixel_color[1u]);
			render_buffer[4u*(w*idy+idx)+2u] = (f32)(pixel_color[2u]);
			render_buffer[4u*(w*idy+idx)+3u] = (f32)(            1.0);
		}
	}
	#endif

	#if 0
	grid_t grid = sim.grid;

	unsigned int grid_dim[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

	unsigned int grid_size = grid_dim[2] * grid_dim[1] * grid_dim[0];

	real bounds_min[3] = {grid.bounds_min[0], grid.bounds_min[1], grid.bounds_min[2]};
	real bounds_max[3] = {grid.bounds_max[0], grid.bounds_max[1], grid.bounds_max[2]};

	real grid_bounds_wid[3] = {
		bounds_max[0] - bounds_min[0],
		bounds_max[1] - bounds_min[1],
		bounds_max[2] - bounds_min[2]
	};

	real* grid_potential = grid_samp(&grid);

	for(unsigned int idx = 0u; idx < w; idx++) {
		for(unsigned int idy = 0u; idy < h; idy++) {
			real value = (real)0.0;

			/*
			for(unsigned int i = 0u; i < grid_dim[2]; i++) {
				int coord[3] = {(int)(idx / 4u), (int)(idy / 4u), (int)i};

				if(0 <= coord[0] && coord[0] < grid_dim[0]
				&& 0 <= coord[1] && coord[1] < grid_dim[1]
				&& 0 <= coord[2] && coord[2] < grid_dim[2]) {
					value += grid_potential[grid_dim[1]*grid_dim[0]*coord[2]+grid_dim[0]*coord[1]+coord[0]];
				}
			}
			*/

			{
				//int coord[3] = {(int)idx, (int)idy, (int)(grid.dim[2] / 2u)};
				int coord[3] = {(int)idx, (int)(grid.dim[1] / 2u), (int)idy};

				if(0 <= coord[0] && coord[0] < grid_dim[0]
				&& 0 <= coord[1] && coord[1] < grid_dim[1]
				&& 0 <= coord[2] && coord[2] < grid_dim[2]) {
					value += grid_potential[grid_dim[1]*grid_dim[0]*coord[2]+grid_dim[0]*coord[1]+coord[0]];
				}
			}

			//value *= -1.000e-14;
			//value *= -1.000e-12;

			real palette_phase = (real)(-1.000e-11) * value;

			real color[3] = {
				(real)0.5*real_cos( (real)TAU * palette_phase - (real)( TAU * (0.0/3.0) ) )+(real)0.5,
				(real)0.5*real_cos( (real)TAU * palette_phase - (real)( TAU * (1.0/3.0) ) )+(real)0.5,
				(real)0.5*real_cos( (real)TAU * palette_phase - (real)( TAU * (2.0/3.0) ) )+(real)0.5
			};

			for(int i = 0; i < 3; i++) {
				color[i] = (real)1.000;
			}

			for(int i = 0; i < 3; i++) {
				color[i] *= (real)(-1.000e-13) * value;
			}

			for(int i = 0; i < 3; i++) {
				color[i] = (real)0.0 <= color[i] ? color[i] : (real)0.0;
				color[i] = (real)1.0 >= color[i] ? color[i] : (real)1.0;
			}

			render_buffer[4u*(w*idy+idx)+0u] = (f32)color[0];
			render_buffer[4u*(w*idy+idx)+1u] = (f32)color[1];
			render_buffer[4u*(w*idy+idx)+2u] = (f32)color[2];
			render_buffer[4u*(w*idy+idx)+3u] = (f32)1.000;
		}
	}
	#endif

	#ifdef LOG_TIMINGS_VIS_DRAW
	LOG_TIMING(log_timings_vis_draw);
	fprintf(log_timings_vis_draw.file, "\n");
	log_sync(&log_timings_vis_draw);
	#endif

	return STARFLOOD_SUCCESS;
}
