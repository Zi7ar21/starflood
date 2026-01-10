// Needed for pthread.h and timing.h
#define _POSIX_C_SOURCE 200809L

#include "visualization.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "config.h"
#include "timing.h"
#include "types.h"

#ifdef ENABLE_VIS_IO_THREAD
#include <pthread.h>
pthread_t vis_io_thread;
pthread_mutex_t vis_io_mutex;
#endif

#ifdef ENABLE_STB_IMAGE_WRITE
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#endif

int image_write_pfm(const void* restrict data, size_t w, size_t h, const char* restrict filename) {
	if(NULL == data || NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	FILE* file = fopen(filename, "wb");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "image_write_pfm()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	// PFM graphic image file format
	// https://netpbm.sourceforge.net/doc/pfm.html
	fprintf(file, "PF\n%zu %zu\n-1.0\n", w, h);

	fwrite(data, sizeof(f32), (size_t)3u * w * h, file);

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose() ", "image_write_pfm()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}

int image_write_ppm(const void* restrict data, size_t w, size_t h, const char* restrict filename) {
	if(NULL == data || NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	FILE* file = fopen(filename, "wb");

	if(NULL == (void*)file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "image_write_ppm()", filename, "wb");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	// PPM Netpbm color image format
	// https://netpbm.sourceforge.net/doc/ppm.html
	fprintf(file, "P6\n%zu %zu %u\n", w, h, 255u);

	fwrite(data, sizeof(unsigned char), (size_t)3u * w * h, file);

	if( 0 != fclose(file) ) {
		fprintf(stderr, "%s error: fclose() ", "image_write_ppm()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}

#ifdef ENABLE_STB_IMAGE_WRITE
int image_write_png(const void* restrict data, size_t w, size_t h, const char* restrict filename) {
	if(NULL == data || NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	stbi_flip_vertically_on_write(0);

	#ifdef IO_PNG_COMPRESSION_LEVEL
	// the default value is 8
	// stbi_zlib_compress() forces quality to be a minimum of 5
	stbi_write_png_compression_level = IO_PNG_COMPRESSION_LEVEL;
	#endif

	int stride = (int)(sizeof(unsigned char) * (size_t)3u * w);

	if( 0 == stbi_write_png(filename, (int)w, (int)h, 3, data, stride) ) {
		fprintf(stderr, "%s error: stbi_write_png() ", "image_write()");
		perror("failed");
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}
#else
int image_write_png(const void* restrict data, size_t w, size_t h, const char* restrict filename) {
	fprintf(stderr, "%s error: Starflood was compiled without PNG image I/O support (in config.h, try uncommenting #define ENABLE_STB_IMAGE_WRITE).\n", "image_write_png()");
	return STARFLOOD_FAILURE;
}
#endif

typedef struct {
	image_filetype_t filetype;

	char filename[STARFLOOD_FILENAME_MAX];

	size_t image_w, image_h; // image dimensions

	void* image_data;
} image_write_arg_t;

int image_write(image_write_arg_t arg) {
	switch(arg.filetype) {
		case IMAGE_FILETYPE_PFM:
			return image_write_pfm(arg.image_data, arg.image_w, arg.image_h, arg.filename);
		case IMAGE_FILETYPE_PPM:
			return image_write_ppm(arg.image_data, arg.image_w, arg.image_h, arg.filename);
		case IMAGE_FILETYPE_PNG:
			return image_write_png(arg.image_data, arg.image_w, arg.image_h, arg.filename);
		default:
			break;
	}

	return STARFLOOD_FAILURE;
}

#ifdef ENABLE_VIS_IO_THREAD
void* image_write_worker(void* arg) {
	if(NULL == arg) {
		return NULL;
	}

	if( 0 != pthread_mutex_lock(&vis_io_mutex) ) {
		fprintf(stderr, "%s error: pthread_mutex_lock(&vis_io_mutex) ", "image_write_worker()");
		perror("failed");
	}

	image_write_arg_t image_write_arg = *(image_write_arg_t*)arg;

	image_write(image_write_arg);

	if( 0 != pthread_mutex_unlock(&vis_io_mutex) ) {
		fprintf(stderr, "%s error: pthread_mutex_unlock(&vis_io_mutex) ", "image_write_worker()");
		perror("failed");
	}

	free(arg);

	return NULL;
}
#endif

int vis_save(const vis_t* vis_ptr, const char* restrict filename, image_filetype_t filetype) {
	TIMING_INIT();

	#ifdef ENABLE_VIS_IO_THREAD
	TIMING_START();

	if( 0 != pthread_join(vis_io_thread, NULL) ) {
	}

	TIMING_STOP();
	TIMING_PRINT("vis_save()", "pthread_join()");
	TIMING_START();

	if( 0 != pthread_mutex_lock(&vis_io_mutex) ) {
	}

	TIMING_STOP();
	TIMING_PRINT("vis_save()", "pthread_mutex_lock()");
	#endif

	vis_t vis = *vis_ptr;

	unsigned int image_w = vis.w;
	unsigned int image_h = vis.h;

	f32* render_buffer = vis.render_buffer;

	TIMING_START();

	if(IMAGE_FILETYPE_PFM == filetype) {
		f32* binary_buffer = vis.binary_buffer;

		// Copy render_buffer->binary_buffer without alpha channel
		for(unsigned int y = 0u; y < image_h; y++) {
			for(unsigned int x = 0u; x < image_w; x++) {
				binary_buffer[3u*(image_w*y+x)+0u] = render_buffer[4u*(image_w*y+x)+0u];
				binary_buffer[3u*(image_w*y+x)+1u] = render_buffer[4u*(image_w*y+x)+1u];
				binary_buffer[3u*(image_w*y+x)+2u] = render_buffer[4u*(image_w*y+x)+2u];
			}
		}
	}

	if(IMAGE_FILETYPE_PPM == filetype || IMAGE_FILETYPE_PNG == filetype) {
		unsigned char* binary_buffer = vis.binary_buffer;

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

				// Gamma 2.2  non-linear transfer function
				//for(int i = 0; i < 3; i++) {
				//	color[i] = (f32)powf( (float)color[i], (float)(1.0/2.2) );
				//}

				if(IMAGE_FILETYPE_PPM == filetype) {
					// BT.709 non-linear transfer function
					for(int i = 0; i < 3; i++) {
						color[i] = (f32)0.018053968510807 <= color[i] ? (f32)1.099296826809442 * (f32)powf( (float)color[i], (float)0.45 ) - (f32)0.099296826809442 : (f32)4.500 * color[i];
					}
				}

				if(IMAGE_FILETYPE_PNG == filetype) {
					// IEC 61966-2-1 sRGB non-linear transfer function
					for(int i = 0; i < 3; i++) {
						color[i] = (f32)0.0031308 < color[i] ? (f32)1.055 * (f32)powf(  (float)color[i], (float)(1.0/2.4) ) - (f32)0.055 : (f32)12.92 * color[i];
					}
				}

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
	}

	TIMING_STOP();
	TIMING_PRINT("vis_save()", "binary_buffer");
	TIMING_START();

	#ifdef ENABLE_VIS_IO_THREAD
	image_write_arg_t* image_write_arg = (image_write_arg_t*)malloc( sizeof(image_write_arg_t) );

	image_write_arg->filetype = filetype;

	strcpy(image_write_arg->filename, filename);

	image_write_arg->image_w = (size_t)image_w;
	image_write_arg->image_h = (size_t)image_h;

	image_write_arg->image_data = vis.binary_buffer;
	#else
	image_write_arg_t image_write_arg;

	image_write_arg.filetype = filetype;

	strcpy(image_write_arg.filename, filename);

	image_write_arg.image_w = (size_t)image_w;
	image_write_arg.image_h = (size_t)image_h;

	image_write_arg.image_data = vis.binary_buffer;
	#endif

	#ifdef ENABLE_VIS_IO_THREAD
	if( 0 != pthread_create(&vis_io_thread, NULL, image_write_worker, (void*)image_write_arg) ) {
		fprintf(stderr, "%s error: pthread_create(&vis_io_thread, NULL, image_write, image_write_arg) ", "vis_save()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("vis_save()", "pthread_create()");
	TIMING_START();

	if( 0 != pthread_mutex_unlock(&vis_io_mutex) ) {
		fprintf(stderr, "%s error: pthread_mutex_unlock(&vis_io_mutex) ", "vis_save()");
		perror("failed");
	}

	TIMING_STOP();
	TIMING_PRINT("vis_save()", "pthread_mutex_unlock()");
	#else
	if( STARFLOOD_SUCCESS != image_write(image_write_arg) ) {
		fprintf(stderr, "%s error: image_write() failed.\n", "vis_save()");
	}

	TIMING_STOP();
	TIMING_PRINT("vis_save()", "image_write()");
	#endif

	return STARFLOOD_SUCCESS;
}
