#pragma once

#include <cstdint>

typedef struct {
	uint32_t z0;
	uint32_t z1;
	uint32_t z2;
	uint32_t z3;
} rng_state_t;

void init_rng(rng_state_t* rng_state);

// floating-point value just less than 2^-32, ensuring the rng never returns >=1.0
const float inv_rng32_max = 2.32830629776081821092e-10f;

void update_rng(rng_state_t* rng_state);
