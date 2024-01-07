#include "rng.h"

#include <stdint.h>

//https://nullprogram.com/blog/2018/07/31/
// exact bias: 0.020888578919738908
uint32_t triple32(uint32_t x) {
	x ^= x >> (uint32_t)17u;
	x *= (uint32_t)0xED5AD4BBu;
	x ^= x >> (uint32_t)11u;
	x *= (uint32_t)0xAC4C1B51u;
	x ^= x >> (uint32_t)15u;
	x *= (uint32_t)0x31848BABu;
	x ^= x >> (uint32_t)14u;
	return x;
}

void update_rng(rng_state_t* rng_state) {
	rng_state->ns = triple32(rng_state->ns);
}

float urand(rng_state_t* rng_state) {
	update_rng(rng_state);

	return 2.32830629776081821092e-10f * (float)rng_state->ns;
}
