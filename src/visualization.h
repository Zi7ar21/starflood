#pragma once

#include <stddef.h>

#include "common.h"
#include "config.h"
#include "simulation.h"
#include "types.h"

// Enumeration of image file types
typedef enum {
	IMAGE_FILETYPE_PFM,
	IMAGE_FILETYPE_PPM,
	IMAGE_FILETYPE_PNG,
	NUM_IMAGE_FILETYPE
} image_filetype_t;

typedef struct {
	unsigned int w;
	unsigned int h;

	void* mem; // Shared memory for the visualization

	i32* atomic_buffer; // atomic buffer for accumulation rasterization
	f32* render_buffer; // render buffer for the post-processed render

	void* binary_buffer; // binary buffer for image file I/O
} vis_t;

int visualization_init(vis_t* restrict visualization, unsigned int w, unsigned int h);

int visualization_free(vis_t* restrict visualization);

int visualization_save(const vis_t* restrict visualization, const char* restrict filename, image_filetype_t filetype);

int visualization_draw(const vis_t* restrict visualization, const sim_t* restrict sim_ptr);
