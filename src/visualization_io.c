// Needed for pthread.h and timing.h
#define _POSIX_C_SOURCE 200809L

#include "visualization.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "common.h"
#include "config.h"
#include "timing.h"
#include "types.h"

#ifdef VISUALIZATION_THREADED_IO
#include <pthread.h>
extern pthread_t vis_io_thread;
extern pthread_mutex_t vis_io_mutex;
#endif

extern struct image_write_param image_write_params;

#if (2 == VISUALIZATION_IMAGE_FORMAT)
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#endif

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
		fprintf(stderr, "%s error: fclose() ", "image_write()");
		perror("failed");
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
	// Copy render_buffer->binary_buffer without alpha channel
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
			unsigned char pixel[3];

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
				color[i] = (f32)0.018053968510807 <= color[i] ? (f32)1.099296826809442 * (f32)powf( (float)color[i], (float)0.45 ) - (f32)0.099296826809442 : (f32)4.500 * color[i];
			}
			#endif

			#if (2 == VISUALIZATION_IMAGE_FORMAT)
			// IEC 61966-2-1 sRGB non-linear transfer function
			for(int i = 0; i < 3; i++) {
				color[i] = (f32)0.0031308 < color[i] ? (f32)1.055 * (f32)powf(  (float)color[i], (float)(1.0/2.4) ) - (f32)0.055 : (f32)12.92 * color[i];
			}
			#endif

			// 8-bit quantization [0, 255]
			{
				for(int i = 0; i < 3; i++) {
					color[i] *= (f32)255.0;
				}

				for(int i = 0; i < 3; i++) {
					color[i] = (f32)roundf( (float)color[i] );
				}

				for(int i = 0; i < 3; i++) {
					//color[i] = (f32)fminf(fmaxf(color[i], (f32)0.0), (f32)1.0);
					color[i] = (f32)(  0.0) < color[i] ? color[i] : (f32)(  0.0);
					color[i] = (f32)(255.0) > color[i] ? color[i] : (f32)(255.0);
				}

				for(int i = 0; i < 3; i++) {
					pixel[i] = (unsigned char)color[i];
				}
			}

			binary_buffer[3u*(image_w*((image_h-1u)-y)+x)+0u] = pixel[0u];
			binary_buffer[3u*(image_w*((image_h-1u)-y)+x)+1u] = pixel[1u];
			binary_buffer[3u*(image_w*((image_h-1u)-y)+x)+2u] = pixel[2u];
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

	#if (2 == VISUALIZATION_IMAGE_FORMAT)
	// the default value is 8
	// stbi_zlib_compress() forces quality to be a minimum of 5
	stbi_write_png_compression_level = 1;

	stbi_flip_vertically_on_write(0);
	#endif

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
