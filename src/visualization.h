#pragma once

#include <stddef.h>

#include "simulation.h"
#include "types.h"

typedef struct {
	unsigned int w;
	unsigned int h;

	void* mem; // Shared memory for the visualization

	i32* atomic_buffer; // atomic buffer for accumulation rasterization
	f32* render_buffer; // render buffer for the post-processed render
} visualization_t;

int visualization_init(visualization_t* visualization, unsigned int w, unsigned int h);

int visualization_draw(const visualization_t* restrict visualization, const simulation_t* restrict simulation);

int visualization_save(const visualization_t* restrict visualization, const char* restrict filename); 

int visualization_free(visualization_t* visualization);
