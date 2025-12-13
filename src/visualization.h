#pragma once

#include <stddef.h>

#include "simulation.h"
#include "types.h"

typedef struct {
	unsigned int w;
	unsigned int h;

	void* mem;

	i32* atomic_buffer;
	f32* render_buffer;
} visualization_t;

int visualization_init(visualization_t* visualization, unsigned int w, unsigned int h);

int visualization_draw(visualization_t* visualization, simulation_t* simulation);

int visualization_save(visualization_t* visualization, const char* restrict filename); 

int visualization_free(visualization_t* visualization);
