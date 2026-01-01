#pragma once

#include <stddef.h>

#include "common.h"
#include "config.h"
#include "simulation.h"
#include "types.h"

typedef struct {
	unsigned int w;
	unsigned int h;

	void* mem; // Shared memory for the visualization

	i32* atomic_buffer; // atomic buffer for accumulation rasterization
	f32* render_buffer; // render buffer for the post-processed render

	// binary image data for file I/O
	#if (0 >= VIS_FILE_FORMAT)
	f32* binary_buffer;
	#else
	unsigned char* binary_buffer;
	#endif
} vis_t;

int visualization_init(vis_t* restrict visualization, unsigned int w, unsigned int h);

int visualization_free(vis_t* restrict visualization);

int visualization_save(const vis_t* restrict visualization, const char* restrict filename);

int visualization_draw(const vis_t* restrict visualization, const sim_t* restrict simulation);

struct image_write_param {
	unsigned int image_w;
	unsigned int image_h;

	#if (0 >= VIS_FILE_FORMAT)
	volatile f32* binary_buffer;
	#else
	volatile unsigned char* binary_buffer;
	#endif

	char filename[STARFLOOD_FILENAME_MAX];
};
