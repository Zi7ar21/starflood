// Needed for posix_memalign()
#define _POSIX_C_SOURCE 200112L

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
#include "types.h"
#include "timing.h"

#ifdef VISUALIZATION_THREADED_IO
#include <pthread.h>
pthread_t vis_io_thread;
pthread_mutex_t vis_io_mutex;
#endif

#if (2 == VISUALIZATION_IMAGE_FORMAT)
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

void* dummy_function(void* arg) {
	return NULL;
}
#endif

struct image_write_param {
	unsigned int image_w;
	unsigned int image_h;

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	volatile f32* binary_buffer;
	#else
	volatile unsigned char* binary_buffer;
	#endif

	char filename[STARFLOOD_FILENAME_MAX];
};

struct image_write_param image_write_params;

#ifdef LOG_TIMINGS_VIS_DRAW
log_t log_timings_vis_draw;
#endif

int visualization_init(vis_t* restrict visualization, unsigned int w, unsigned int h) {
	#ifdef VISUALIZATION_THREADED_IO
	if( 0 != pthread_create(&vis_io_thread, NULL, dummy_function, NULL) ) {
		fprintf(stderr, "%s error: pthread_create(&vis_io_thread, NULL, dummy_function, NULL) ", "visualization_init()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	if( 0 != pthread_mutex_init(&vis_io_mutex, NULL) ) {
		fprintf(stderr, "%s error: pthread_create(&vis_io_thread, NULL, dummy_function, NULL) ", "visualization_init()");
		perror("failed");
		pthread_join(vis_io_thread, NULL);
		return STARFLOOD_FAILURE;
	}
	#endif

	TIMING_INIT();

	#ifdef LOG_TIMINGS_VIS_DRAW
	if( STARFLOOD_SUCCESS != log_init(&log_timings_vis_draw, LOG_TIMINGS_VIS_DRAW) ) {
		return STARFLOOD_FAILURE;
	}

	fprintf(log_timings_vis_draw.file, "%s,%s,%s,%s,%s\n", "step_number", "clear_atomic_buffer", "calc_matMVP", "rasterize_atomic", "finalize");

	fflush(log_timings_vis_draw.file);
	#endif

	vis_t vis = *visualization;

	void* mem = NULL;

	i32* atomic_buffer = (i32*)NULL;
	f32* render_buffer = (f32*)NULL;

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	f32* binary_buffer = (f32*)NULL;
	#else
	unsigned char* binary_buffer = (unsigned char*)NULL;
	#endif

	if( sizeof(i32) != sizeof(f32) ) {
		return STARFLOOD_FAILURE;
	}

	printf("Visualization Memory Addresses:\n");

	size_t atomic_buffer_length = (size_t)4u * (size_t)w * (size_t)h;
	size_t render_buffer_length = (size_t)4u * (size_t)w * (size_t)h;
	size_t binary_buffer_length = (size_t)3u * (size_t)w * (size_t)h;

	size_t atomic_buffer_offset = (size_t)0u;
	size_t render_buffer_offset = atomic_buffer_offset + atomic_buffer_length;
	size_t binary_buffer_offset = render_buffer_offset + render_buffer_length;

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
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

	printf("  mem: %p\n", mem);

	TIMING_START();

	memset(mem, 0, mem_size);

	TIMING_STOP();
	TIMING_PRINT("visualization_init()", "memset()");

	atomic_buffer = (i32*)mem + (size_t)atomic_buffer_offset;
	render_buffer = (f32*)mem + (size_t)render_buffer_offset;

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	binary_buffer = (f32*)mem + (size_t)binary_buffer_offset;
	#else
	binary_buffer = (unsigned char*)mem + (size_t)binary_buffer_offset;
	#endif

	printf("  atomic_buffer: %p (+%zu)\n", (void*)atomic_buffer, atomic_buffer_offset);
	printf("  render_buffer: %p (+%zu)\n", (void*)render_buffer, render_buffer_offset);
	printf("  binary_buffer: %p (+%zu)\n", (void*)binary_buffer, binary_buffer_offset);
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
		#if (0 >= VISUALIZATION_IMAGE_FORMAT)
		binary_buffer[i] = (f32)0.000;
		#else
		binary_buffer[i] = (unsigned char)0u;
		#endif
	}
	TIMING_STOP();
	TIMING_PRINT("visualization_init()", "initialize binary_buffer");

	#if (2 == VISUALIZATION_IMAGE_FORMAT)
	stbi_flip_vertically_on_write(0);
	stbi_write_png_compression_level = 6;
	#endif

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
	#ifdef VISUALIZATION_THREADED_IO
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

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	vis.binary_buffer = (f32*)NULL;
	#else
	vis.binary_buffer = (unsigned char*)NULL;
	#endif

	*visualization = vis;

	return STARFLOOD_SUCCESS;
}

void* image_write(void* arg) {
	#ifdef VISUALIZATION_THREADED_IO
	pthread_mutex_lock(&vis_io_mutex);
	#endif

	unsigned int image_w = image_write_params.image_w;
	unsigned int image_h = image_write_params.image_h;

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	f32* binary_buffer = (f32*)image_write_params.binary_buffer;
	#else
	unsigned char* binary_buffer = (unsigned char*)image_write_params.binary_buffer;
	#endif

	char* filename = image_write_params.filename;

	#if (2 == VISUALIZATION_IMAGE_FORMAT)
	int stride = (int)(sizeof(unsigned char) * (size_t)3u * (size_t)image_w);

	if( 0 == stbi_write_png(filename, (int)image_w, (int)image_h, 3, binary_buffer, stride) ) {
	}
	#endif

	#if (1 >= VISUALIZATION_IMAGE_FORMAT)
	FILE* file = fopen(filename, "wb");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "image_write()", filename, "wb");
		perror("failed");
		#ifdef VISUALIZATION_THREADED_IO
		pthread_mutex_unlock(&vis_io_mutex);
		#endif
		return NULL;
	}

	#if (0 == VISUALIZATION_IMAGE_FORMAT)
	// PFM graphic image file format
	// https://netpbm.sourceforge.net/doc/pfm.html
	fprintf(file, "PF\n%zu %zu\n-1.0\n", (size_t)image_w, (size_t)image_h);

	fwrite(binary_buffer, sizeof(f32), (size_t)3u * (size_t)image_w * (size_t)image_h, file);
	#endif

	#if (1 == VISUALIZATION_IMAGE_FORMAT)
	// PPM Netpbm color image format
	// https://netpbm.sourceforge.net/doc/ppm.html
	fprintf(file, "P6\n%zu %zu %u\n", (size_t)image_w, (size_t)image_h, 255u);

	fwrite(binary_buffer, sizeof(unsigned char), (size_t)3u * (size_t)image_w * (size_t)image_h, file);
	#endif

	if( 0 != fclose(file) ) {
		#ifdef VISUALIZATION_THREADED_IO
		pthread_mutex_unlock(&vis_io_mutex);
		#endif
		return NULL;
		//return STARFLOOD_FAILURE;
	}
	#endif

	#ifdef VISUALIZATION_THREADED_IO
	pthread_mutex_unlock(&vis_io_mutex);
	#endif

	return NULL;
}

int visualization_save(const vis_t* restrict visualization, const char* restrict filename) {
	#ifdef VISUALIZATION_THREADED_IO
	pthread_mutex_lock(&vis_io_mutex);
	pthread_join(vis_io_thread, NULL);
	#endif

	TIMING_INIT();

	vis_t vis = *visualization;

	unsigned int image_w = vis.w;
	unsigned int image_h = vis.h;

	f32* render_buffer = vis.render_buffer;

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	f32* binary_buffer = vis.binary_buffer;
	#else
	unsigned char* binary_buffer = vis.binary_buffer;
	#endif

	TIMING_START();

	#if (0 >= VISUALIZATION_IMAGE_FORMAT)
	for(unsigned int y = 0u; y < image_h; y++) {
		for(unsigned int x = 0u; x < image_w; x++) {
			binary_buffer[3u*(image_w*y+x)+0u] = render_buffer[4u*(image_w*y+x)+0u];
			binary_buffer[3u*(image_w*y+x)+1u] = render_buffer[4u*(image_w*y+x)+1u];
			binary_buffer[3u*(image_w*y+x)+2u] = render_buffer[4u*(image_w*y+x)+2u];
		}
	}
	#else
	for(unsigned int y = 0u; y < image_h; y++) {
		for(unsigned int x = 0u; x < image_w; x++) {
			f32 color[3] = {
				(f32)render_buffer[4u*(image_w*y+x)+0u],
				(f32)render_buffer[4u*(image_w*y+x)+1u],
				(f32)render_buffer[4u*(image_w*y+x)+2u]
			};

			// Clamp values to range [0, 1]
			for(int i = 0; i < 3; i++) {
				//color[i] = (f32)fminf(fmaxf(color[i], (f32)0.0), (f32)1.0);

				color[i] = (f32)0.0 < color[i] ? color[i] : (f32)0.0;
				color[i] = (f32)1.0 > color[i] ? color[i] : (f32)1.0;
			}

			/*
			// Gamma 2.2  non-linear transfer function
			for(int i = 0; i < 3; i++) {
				color[i] = (f32)powf( color[i], (f32)(1.0/2.2) );
			}
			*/

			#if (1 == VISUALIZATION_IMAGE_FORMAT)
			// BT.709 non-linear transfer function
			for(int i = 0; i < 3; i++) {
				color[i] = (f32)0.018053968510807 <= color[i] ? (f32)1.099296826809442 * (f32)powf(color[i], (f32)0.45) - (f32)0.099296826809442 : (f32)4.500 * color[i];
			}
			#endif

			#if (2 == VISUALIZATION_IMAGE_FORMAT)
			// IEC 61966-2-1 sRGB non-linear transfer function
			for(int i = 0; i < 3; i++) {
				color[i] = (f32)0.0031308 < color[i] ? (f32)1.055 * (f32)powf( color[i], (f32)(1.0/2.4) ) - (f32)0.055 : (f32)12.92 * color[i];
			}
			#endif

			// 8-bit quantization [0, 255]
			unsigned char rgb[3] = {
				(unsigned char)( (f32)255.0 * color[0] ),
				(unsigned char)( (f32)255.0 * color[1] ),
				(unsigned char)( (f32)255.0 * color[2] )
			};

			binary_buffer[3u*(image_w*((image_h-1u)-y)+x)+0u] = rgb[0u];
			binary_buffer[3u*(image_w*((image_h-1u)-y)+x)+1u] = rgb[1u];
			binary_buffer[3u*(image_w*((image_h-1u)-y)+x)+2u] = rgb[2u];
		}
	}
	#endif

	TIMING_STOP();
	TIMING_PRINT("visualization_save()", "binary_buffer");
	TIMING_START();

	image_write_params.image_w = image_w;
	image_write_params.image_h = image_h;

	image_write_params.binary_buffer = (void*)binary_buffer;

	strcpy(image_write_params.filename, filename);

	#ifdef VISUALIZATION_THREADED_IO
	pthread_create(&vis_io_thread, NULL, image_write, NULL);
	#else
	image_write(NULL);
	#endif

	TIMING_STOP();
	TIMING_PRINT("visualization_save()", "image_write()");

	#ifdef VISUALIZATION_THREADED_IO
	pthread_mutex_unlock(&vis_io_mutex);
	#endif

	return STARFLOOD_SUCCESS;
}

int visualization_draw(const vis_t* restrict visualization, const sim_t* restrict simulation) {
	TIMING_INIT();

	sim_t sim = *simulation;

	vis_t vis = *visualization;

	unsigned int w = vis.w;
	unsigned int h = vis.h;

	i32* atomic_buffer = vis.atomic_buffer;
	f32* render_buffer = vis.render_buffer;

	unsigned int step_number = sim.step_number;

	unsigned int N = sim.N;

	real* pos = sim.pos;
	real* vel = sim.vel; // for motion blur

	real* pot = sim.pot; // for color

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

		//bet = (real)0.0000 * (real)TAU;
		//bet = (real)0.0625 * (real)TAU;
		//bet = (real)0.1250 * (real)TAU;
		//bet = (real)0.2500 * (real)TAU;
		bet = (real)1.000e-4 * time * (real)TAU;

		//gam = (real)0.0000 * (real)TAU;
		//gam = (real)0.0625 * (real)TAU;
		//gam = (real)0.1250 * (real)TAU;
		//gam = (real)0.2500 * (real)TAU;

		#ifdef STARFLOOD_DOUBLE_PRECISION
		real cos_alp = (real)cos(alp), sin_alp = (real)sin(alp);
		real cos_bet = (real)cos(bet), sin_bet = (real)sin(bet);
		real cos_gam = (real)cos(gam), sin_gam = (real)sin(gam);
		#else
		real cos_alp = (real)cosf(alp), sin_alp = (real)sinf(alp);
		real cos_bet = (real)cosf(bet), sin_bet = (real)sinf(bet);
		real cos_gam = (real)cosf(gam), sin_gam = (real)sinf(gam);
		#endif

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
	double ortho_scale = exp2(ORTHO_SCALE);
	#else
	double ortho_scale = 1.0;
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
		#ifdef ENABLE_OFFLOADING
		#pragma omp target teams distribute parallel for map(tofrom: atomic_buffer[:4u*w*h]) map(to: pos[:3u*N], vel[:3u*N], pot[:N])
		#else
		#pragma omp parallel for schedule(dynamic, 1024)
		#endif
	#endif
	for(unsigned int idx = 0u; idx < N; idx++) {
		// https://cs418.cs.illinois.edu/website/text/math2.html

		// homogenous coordinates
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

		real color_phase = (real)3.000e0 * pot[idx];

		real color[3] = {
			#ifdef STARFLOOD_DOUBLE_PRECISION
			(real)0.5*(real)cos( (real)TAU * color_phase - (real)( (0.0/3.0) * TAU ) )+(real)0.5,
			(real)0.5*(real)cos( (real)TAU * color_phase - (real)( (1.0/3.0) * TAU ) )+(real)0.5,
			(real)0.5*(real)cos( (real)TAU * color_phase - (real)( (2.0/3.0) * TAU ) )+(real)0.5
			#else
			(real)0.5*(real)cosf( (real)TAU * color_phase - (real)( (0.0/3.0) * TAU ) )+(real)0.5,
			(real)0.5*(real)cosf( (real)TAU * color_phase - (real)( (1.0/3.0) * TAU ) )+(real)0.5,
			(real)0.5*(real)cosf( (real)TAU * color_phase - (real)( (2.0/3.0) * TAU ) )+(real)0.5
			#endif
		};

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

		for(unsigned int samp_number = 0u; samp_number < (unsigned int)SPATIAL_SAMPLES; samp_number++) {
			real p_i[3] = {
				p[0],
				p[1],
				p[2]
			};

			u32 s[4] = {(u32)samp_number, (u32)idx, (u32)step_number, (u32)0xCF52CA01u};

			pcg4d(s);
			//pcg4d(s); // second round for better statistical quality

			real time_offset = (real)0.0;

			real sample_offset[2] = {(real)0.0, (real)0.0};

			i32 accumulation_color[4] = {
				(i32)( (real)1.000e2 * color[0] ),
				(i32)( (real)1.000e2 * color[1] ),
				(i32)( (real)1.000e2 * color[2] ),
				(i32)1
			};

			{
				pcg4d(s);

				real r[4] = {
					(real)INV_PCG32_MAX * (real)s[0],
					(real)INV_PCG32_MAX * (real)s[1],
					(real)INV_PCG32_MAX * (real)s[2],
					(real)INV_PCG32_MAX * (real)s[3]
				};

				// Box-Muller Transform
				// https://en.wikipedia.org/wiki/Box–Muller_transform
				real n[2] = {
					#ifdef STARFLOOD_DOUBLE_PRECISION
					(real)sqrt( (real)(-2.0) * (real)log(r[2]) ) * (real)cos( (real)TAU * r[3] ),
					(real)sqrt( (real)(-2.0) * (real)log(r[2]) ) * (real)sin( (real)TAU * r[3] )
					#else
					(real)sqrtf( (real)(-2.0) * (real)logf(r[2]) ) * (real)cosf( (real)TAU * r[3] ),
					(real)sqrtf( (real)(-2.0) * (real)logf(r[2]) ) * (real)sinf( (real)TAU * r[3] )
					#endif
				};

				// spatial anti-aliasing (gaussian)
				{
					const real filter_sigma = (real)0.800;

					sample_offset[0] = filter_sigma * n[0];
					sample_offset[1] = filter_sigma * n[1];
				}

				{
					#ifdef SHUTTER_SPEED
					time_offset += (real)SHUTTER_SPEED * (r[1] - (real)0.5); // motion blur
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
				double x[2], y[2];

				x[0] = p_i[1]; // Y
				x[1] = p_i[2]; // Z

				//const double theta = 0.0000 * TAU; // side-on
				const double theta = 0.0625 * TAU; // slightly tilted
				//const double theta = 0.1250 * TAU; // edge-on
				//const double theta = 0.2500 * TAU; // top-down

				double cos_theta = cos(theta), sin_theta = sin(theta);

				y[0] = cos_theta*x[0]-sin_theta*x[1];
				y[1] = sin_theta*x[0]+cos_theta*x[1];

				p_i[1] = y[0];
				p_i[2] = y[1];
			}

			// rotate2d
			{
				double x[2], y[2];

				x[0] = p_i[0]; // X
				x[1] = p_i[2]; // Z

				//const double theta = 0.0000 * TAU;
				//const double theta = 0.0625 * TAU;
				//const double theta = 0.1250 * TAU;
				//const double theta = 0.2500 * TAU;
				double theta = 0.001 * (time + time_offset) * TAU;

				double cos_theta = cos(theta), sin_theta = sin(theta);

				y[0] = cos_theta*x[0]-sin_theta*x[1];
				y[1] = sin_theta*x[0]+cos_theta*x[1];

				p_i[0] = y[0];
				p_i[2] = y[1];
			}
			*/

			int coord[2] = {
				(int)( ( (real)h * ( (real)0.5 * p_i[0] ) ) + ( (real)0.5 * (real)w ) + sample_offset[0] - (real)0.5),
				(int)( ( (real)h * ( (real)0.5 * p_i[1] ) ) + ( (real)0.5 * (real)h ) + sample_offset[1] - (real)0.5)
			};

			if( (int)0 <= coord[0] && coord[0] < (int)w && (int)0 <= coord[1] && coord[1] < (int)h ) {
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

	// Since atomic_buffer colors are quantized integers, it needs a floating-point scaling factor
	double pixel_value_scale = 1.0;

	pixel_value_scale *= 1.000e-2;

	#ifdef SPATIAL_SAMPLES
	pixel_value_scale *= 1.0 / (double)SPATIAL_SAMPLES;
	#endif

	#ifdef EXPOSURE
	pixel_value_scale *= exp2(EXPOSURE);
	#endif

	TIMING_START();

	// read the atomic buffer and finish rendering the visualization
	for(unsigned int i = 0u; i < w * h; i++) {
		double rgba[4] = {
			(double)0.250,
			(double)0.250,
			(double)0.250,
			(double)1.000
		};

		for(unsigned int j = 0u; j < 4u; j++) {
			i32 val = (i32)0;

			#ifdef _OPENMP
			#pragma omp atomic read
			#endif
			val = atomic_buffer[4u*i+j];

			rgba[j] = pixel_value_scale * (double)val;
		}

		// tonemapping
		for(unsigned int j = 0u; j < 4u; j++) {
			rgba[j] = tanh(rgba[j]);
		}

		// set alpha channel to 1.000
		rgba[3u] = (double)1.000;

		for(unsigned int j = 0u; j < 4u; j++) {
			render_buffer[4u*i+j] = (f32)rgba[j];
		}
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_draw()", "finalize");
	#ifdef LOG_TIMINGS_VIS_DRAW
	LOG_TIMING(log_timings_vis_draw);
	fprintf(log_timings_vis_draw.file, "\n");
	fflush(log_timings_vis_draw.file);
	#endif

	return STARFLOOD_SUCCESS;
}
