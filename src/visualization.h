#pragma once

#include <stddef.h>

#include "simulation.h"

typedef struct {
	int id;

	unsigned int w;
	unsigned int h;

	int* atomic_buffer;

	float* render_buffer;
} visualization_t;

int visualization_init(visualization_t* visualization, unsigned int w, unsigned int h);

int visualization_draw(visualization_t* visualization, simulation_t* simulation);

int visualization_free(visualization_t* visualization);
