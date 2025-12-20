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
#include "rng.h"
#include "simulation.h"
#include "types.h"
#include "timing.h"

int visualization_init(visualization_t* restrict visualization, unsigned int w, unsigned int h) {
	TIMING_INIT();

	visualization_t vis = *visualization;

	void* mem = NULL;

	i32* atomic_buffer = (i32*)NULL;
	f32* render_buffer = (f32*)NULL;

	if( sizeof(i32) != sizeof(f32) ) {
		return STARFLOOD_FAILURE;
	}

	printf("Visualization Memory Addresses:\n");

	size_t atomic_buffer_length = (size_t)4u * (size_t)w * (size_t)h;
	size_t render_buffer_length = (size_t)4u * (size_t)w * (size_t)h;

	size_t atomic_buffer_offset = (size_t)0u;
	size_t render_buffer_offset = atomic_buffer_offset + atomic_buffer_length;

	size_t mem_size = (sizeof(i32)*atomic_buffer_length)+(sizeof(i32)*render_buffer_length);

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
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

	memset(mem, 0, mem_size);

	atomic_buffer = (i32*)mem + (size_t)atomic_buffer_offset;
	render_buffer = (f32*)mem + (size_t)render_buffer_offset;

	printf("  atomic_buffer: %p (+%zu)\n", (void*)atomic_buffer, atomic_buffer_offset);
	printf("  render_buffer: %p (+%zu)\n", (void*)render_buffer, render_buffer_offset);
	printf("\n");

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		atomic_buffer[i] = (i32)0;
	}

	for(unsigned int i = 0u; i < 4u * w * h; i++) {
		render_buffer[i] = (f32)0.000;
	}

	vis.w = w;
	vis.h = h;
	vis.mem = mem;
	vis.atomic_buffer = atomic_buffer;
	vis.render_buffer = render_buffer;

	*visualization = vis;

	return STARFLOOD_SUCCESS;
}

int visualization_free(visualization_t* restrict visualization) {
	TIMING_INIT();

	visualization_t vis = *visualization;

	free(vis.mem);

	vis.mem = NULL;

	vis.atomic_buffer = (i32*)NULL;
	vis.render_buffer = (f32*)NULL;

	*visualization = vis;

	return STARFLOOD_SUCCESS;
}

int visualization_save(const visualization_t* restrict visualization, const char* restrict filename) {
	TIMING_INIT();

	visualization_t vis = *visualization;

	unsigned int image_w = vis.w;
	unsigned int image_h = vis.h;

	f32* render_buffer = vis.render_buffer;

	TIMING_START();

	FILE* file = fopen(filename, "wb");

	TIMING_STOP();
	TIMING_PRINT("visualization_save()", "fopen()");
	TIMING_START();

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(%s, \"%s\") ", "visualization_save()", filename, "wb");

		perror("failed");

		return STARFLOOD_FAILURE;
	}

	TIMING_START();

	/*
	// PFM graphic image file format
	// https://netpbm.sourceforge.net/doc/pfm.html
	fprintf(file, "PF\n%zu %zu\n-1.0\n", (size_t)image_w, (size_t)image_h);

	for(unsigned int y = 0u; y < image_h; y++) {
		for(unsigned int x = 0u; x < image_w; x++) {
			f32 rgb[3] = {
				render_buffer[4u*(image_w*y+x)+0u],
				render_buffer[4u*(image_w*y+x)+1u],
				render_buffer[4u*(image_w*y+x)+2u]
			};

			fwrite(rgb, sizeof(f32), (size_t)3u, file);
		}
	}
	*/

	// PPM Netpbm color image format
	// https://netpbm.sourceforge.net/doc/ppm.html
	fprintf(file, "P6\n%zu %zu %u\n", (size_t)image_w, (size_t)image_h, 255u);

	for(unsigned int y = 0u; y < image_h; y++) {
		for(unsigned int x = 0u; x < image_w; x++) {
			f32 color[3] = {
				(f32)render_buffer[4u*(image_w*((image_h-1u)-y)+x)+0u],
				(f32)render_buffer[4u*(image_w*((image_h-1u)-y)+x)+1u],
				(f32)render_buffer[4u*(image_w*((image_h-1u)-y)+x)+2u]
			};

			// Clamp values to range [0, 1]
			for(int i = 0; i < 3; i++) {
				color[i] = fminf(fmaxf(color[i], 0.0), 1.0);
			}

			// BT.709 non-linear encoding
			for(int i = 0; i < 3; i++) {
				color[i] = (f32)0.018053968510807 <= color[i] ? (f32)1.099296826809442 * powf(color[i], (f32)0.45) - (f32)0.099296826809442 : (f32)4.500 * color[i];
			}

			// 8-bit quantization [0, 255]
			unsigned char rgb[3] = {
				(unsigned char)( (f32)255.0 * color[0] ),
				(unsigned char)( (f32)255.0 * color[1] ),
				(unsigned char)( (f32)255.0 * color[2] )
			};

			fwrite(rgb, sizeof(unsigned char), (size_t)3u, file);
		}
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_save()", "fwrite()");
	TIMING_START();

	if( 0 != fclose(file) ) {
		return STARFLOOD_FAILURE;
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_save()", "fclose()");

	return STARFLOOD_SUCCESS;
}

int visualization_draw(const visualization_t* restrict visualization, const simulation_t* restrict simulation) {
	TIMING_INIT();

	simulation_t sim = *simulation;

	visualization_t vis = *visualization;

	unsigned int w = vis.w;
	unsigned int h = vis.h;

	i32* atomic_buffer = vis.atomic_buffer;
	f32* render_buffer = vis.render_buffer;

	unsigned int step_number = sim.step_number;

	unsigned int N = sim.N;

	real* pot = sim.pot; // for color
	real* pos = sim.pos;
	real* vel = sim.vel; // for motion blur

	#ifdef OUTPUT_INTERVAL
	double time = TIMESTEP_SIZE * (double)step_number;
	#else
	double time = TIMESTEP_SIZE * (double)step_number;
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
	TIMING_START();

	// model matrix
	double matM[16];

	// model matrix = rotation matrix
	{
		// I feel genuinely sorry for the compiler
		// that has to optimize all of this lol

		double alp = 0.0; // roll
		double bet = 0.0; // pitch
		double gam = 0.0; // yaw

		//bet = 0.0000 * TAU;
		bet = 0.0625 * TAU;
		//bet = 0.1250 * TAU;
		//bet = 0.2500 * TAU;
		//bet = 0.001 * time * TAU;

		//gam = 0.0000 * TAU;
		//gam = 0.0625 * TAU;
		gam = 0.1250 * TAU;
		//gam = 0.2500 * TAU;

		double cos_alp = cos(alp), sin_alp = sin(alp);
		double cos_bet = cos(bet), sin_bet = sin(bet);
		double cos_gam = cos(gam), sin_gam = sin(gam);

		double rot_mat[16] = {
			(cos_alp*cos_bet), (cos_alp*sin_bet*sin_gam-sin_alp*cos_gam), (cos_alp*sin_bet*cos_gam+sin_alp*sin_gam), (0.0),
			(sin_alp*cos_bet), (sin_alp*sin_bet*sin_gam+cos_alp*cos_gam), (sin_alp*sin_bet*cos_gam-cos_alp*sin_gam), (0.0),
			(       -sin_bet), (                        cos_bet*sin_gam), (                        cos_bet*cos_gam), (0.0),
			(            0.0), (                                    0.0), (                                    0.0), (1.0)
		};

		for(int i = 0; i < 16; i++) {
			matM[i] = rot_mat[i];
		}
	}

	// view matrix
	double matV[16] = {
		( 1.000), ( 0.000), ( 0.000), ( 0.000),
		( 0.000), ( 1.000), ( 0.000), ( 0.000),
		( 0.000), ( 0.000), ( 1.000), ( 0.000),
		( 0.000), ( 0.000), ( 0.000), ( 1.000)
	};

	#ifdef ORTHO_SCALE
	double ortho_scale = exp2(ORTHO_SCALE);
	#else
	double ortho_scale = 1.0;
	#endif

	// projection matrix
	double matP[16] = {
		(ortho_scale), (        0.0), (        0.0), (        0.0),
		(        0.0), (ortho_scale), (        0.0), (        0.0),
		(        0.0), (        0.0), (ortho_scale), (        0.0),
		(        0.0), (        0.0), (        0.0), (        1.0)
	};

	double matMVP[16];

	// C = A * B
	{
		double A[16], B[16], C[16];

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
		double A[16], B[16], C[16];

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

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOADING
		#pragma omp target teams distribute parallel for map(tofrom: atomic_buffer[:4u*w*h]) map(to: pot[:N], pos[:3u*N], vel[:3u*N])
		#else
		#pragma omp parallel for schedule(dynamic, 1024)
		#endif
	#endif
	for(unsigned int idx = 0u; idx < N; idx++) {
		// https://cs418.cs.illinois.edu/website/text/math2.html

		// homogenous coordinates
		double p[4] = {
			(double)(pos[3u*idx+0u]),
			(double)(pos[3u*idx+1u]),
			(double)(pos[3u*idx+2u]),
			(double)(           1.0)
		};

		double v[3] = {
			(double)(vel[3u*idx+0u]),
			(double)(vel[3u*idx+1u]),
			(double)(vel[3u*idx+2u])
		};

		double color[3] = {
			0.5 * cos( TAU * ( (-10000.0 * (double)pot[idx])-(0.0/3.0) ) ) + 0.5,
			0.5 * cos( TAU * ( (-10000.0 * (double)pot[idx])-(1.0/3.0) ) ) + 0.5,
			0.5 * cos( TAU * ( (-10000.0 * (double)pot[idx])-(2.0/3.0) ) ) + 0.5,
		};

		// y = A * x
		{
			double A[16], x[4], y[4];

			for(int i = 0; i < 16; i++) {
				A[i] = matMVP[i];
			}

			for(int i = 0; i < 16; i++) {
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
			double p_i[3] = {
				p[0],
				p[1],
				p[2]
			};

			u32 s[4] = {(u32)samp_number, (u32)idx, (u32)step_number, (u32)0xCF52CA01u};

			pcg4d(s);
			//pcg4d(s); // second round for better statistical quality

			double time_offset = 0.0;

			double sample_offset[2] = {0.0, 0.0};

			i32 accumulation_color[4] = {
				(i32)(100.0 * color[0]),
				(i32)(100.0 * color[1]),
				(i32)(100.0 * color[2]),
				(i32)1
			};

			{
				pcg4d(s);

				double r[4] = {
					INV_PCG32_MAX * (double)s[0],
					INV_PCG32_MAX * (double)s[1],
					INV_PCG32_MAX * (double)s[2],
					INV_PCG32_MAX * (double)s[3]
				};

				// spatial anti-aliasing (gaussian)
				// https://en.wikipedia.org/wiki/Box–Muller_transform
				{
					sample_offset[0] = 0.750 * sqrt( -2.0 * log(r[2]) ) * cos(TAU * r[3]);
					sample_offset[1] = 0.750 * sqrt( -2.0 * log(r[2]) ) * sin(TAU * r[3]);
				}

				{
					#ifdef SHUTTER_SPEED
					time_offset += SHUTTER_SPEED * (r[1] - 0.5); // motion blur
					#endif

					time_offset *= (double)OUTPUT_INTERVAL * TIMESTEP_SIZE;

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
				(int)( ( (double)h * (0.500 * p_i[0]) ) + (0.5 * (double)w) + sample_offset[0] - 0.5),
				(int)( ( (double)h * (0.500 * p_i[1]) ) + (0.5 * (double)h) + sample_offset[1] - 0.5)
			};

			if((int)0 <= coord[0] && coord[0] < (int)w && (int)0 <= coord[1] && coord[1] < (int)h) {
				size_t pixel_index = (size_t)w*(size_t)coord[1]+(size_t)coord[0];

				for(size_t i = (size_t)0u; i < (size_t)4u; i++) {
					#ifdef _OPENMP
					#pragma omp atomic update
					#endif
					atomic_buffer[4u*pixel_index+i] += accumulation_color[i];
				}
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("visualization_draw()", "atomic_rasterization");
	TIMING_START();

	// Since atomic_buffer colors are quantized integers, it needs a floating-point scaling factor
	double pixel_value_scale = 1.0;

	pixel_value_scale *= 0.01;

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

	return STARFLOOD_SUCCESS;
}
