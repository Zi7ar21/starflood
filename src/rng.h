#pragma once

#include <stdint.h>

typedef struct {
	uint32_t ns;
} rng_state_t;

void update_rng(rng_state_t* rng_state);

float urand(rng_state_t* rng_state);
