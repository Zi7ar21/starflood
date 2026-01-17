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
	unsigned int sizex, sizey; // Visualization dimensions

	void* mem; // Shared memory for the visualization

	i32* atomic_buffer; // atomic buffer for accumulation rasterization
	f32* render_buffer; // render buffer for the post-processed render

	void* binary_buffer; // binary buffer for image file I/O
} vis_t;

// Initializes a visualization with with dimensions (sizex, sizey)
int vis_init(vis_t* restrict vis_ptr, const unsigned int sizex, const unsigned int sizey);

// Frees resources allocated by a visualization
int vis_free(vis_t* restrict vis_ptr);

// Saves the visualization to a file
int vis_save(const vis_t* restrict vis_ptr, const char* restrict filename, image_filetype_t filetype);

// Update the visualization by drawing a simulation
int vis_draw(const vis_t* restrict vis_ptr, const sim_t* restrict sim_ptr);
