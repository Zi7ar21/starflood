#pragma once

#include <cstdint>

#include "constants.h"
#include "types.hpp"

// https://www.shadertoy.com/view/NldfDn

// ##### Random Number Generator #####

#ifdef USE_PCG_RANDOM
// See michael0884's usage of PCG Random
// https://www.shadertoy.com/view/wltcRS
// https://www.shadertoy.com/view/WttyWX

struct rng_state_t {
	uvec4 seed;
};

void init_rng(inout rng_state_t rng_state, vec2 pixel_coord, vec2 res, int frame_number) {
	uvec2 coord = uvec2(pixel_coord);

	rng_state.seed = uvec4(coord, uint(frame_number), coord.x + coord.y);

	//rng_state.seed = uvec4(frame, frame * 15843, frame * 31 + 4566, frame * 2345 + 58585);
}

// http://www.pcg-random.org/
uvec4 pcg4d(uvec4 v) {
	v = v * 0x0019660Du + 0x3C6EF35Fu;
	v.x += v.y * v.w; v.y += v.z * v.x; v.z += v.x * v.y; v.w += v.y * v.z;
	v = v ^ (v >> 16u);
	v.x += v.y * v.w; v.y += v.z * v.x; v.z += v.x * v.y; v.w += v.y * v.z;

	return v;
}

float rand(inout rng_state_t rng_state) {
	rng_state.seed = pcg4d(rng_state.seed);

	return INV_RNG32_MAX * float(rng_state.seed.x);
}

vec2 rand2(inout rng_state_t rng_state) {
	rng_state.seed = pcg4d(rng_state.seed);

	return INV_RNG32_MAX * vec2(rng_state.seed.xy);
}

vec3 rand3(inout rng_state_t rng_state) {
	rng_state.seed = pcg4d(rng_state.seed);

	return INV_RNG32_MAX * vec3(rng_state.seed.xyz);
}

vec4 rand4(inout rng_state_t rng_state) {
	rng_state.seed = pcg4d(rng_state.seed);

	return INV_RNG32_MAX * vec4(rng_state.seed);
}
#else
typedef struct {
	uint32_t x;
} rng_state_t;

/*
void init_rng(rng_state_t* rng_state, uint32_t coord_x, uint32_t coord_y,vec2 pixel_coord, vec2 res, int frame_number) {
	pixel_coord = floor(pixel_coord);

	rng_state.seed = uint(frame_number) * uint(res.x * res.y) + uint(pixel_coord.x + res.x * pixel_coord.y) + 1u;
}
*/

// Triple32 Hash: https://nullprogram.com/blog/2018/07/31/
// this hash in particular is remarkable because it's a statistically perfect 32-bit integer hash (of this kind)
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

// uniform random
float urand1(rng_state_t* rng_state) {
	rng_state->x = triple32(rng_state->x);
	return INV_RNG32_MAX * (float)rng_state->x;
}

/*
vec2 urand2(inout rng_state_t rng_state) {
	vec2 v;
	v.x = urand1(rng_state);
	v.y = urand1(rng_state);
	return v;
}

vec3 urand3(inout rng_state_t rng_state) {
	vec3 v;
	v.x = rand(rng_state);
	v.y = rand(rng_state);
	v.z = rand(rng_state);
	return v;
}

vec4 urand4(inout rng_state_t rng_state) {
	vec4 v;
	v.x = rand(rng_state);
	v.y = rand(rng_state);
	v.z = rand(rng_state);
	v.w = rand(rng_state);
	return v;
}
*/
#endif

/*
vec2 udir2(inout rng_state_t rng_state) {
	float z = rand(rng_state);
	float t = TAU * z;
	float c = cos(t), s = sin(t);
	return vec2(c, s);
}

vec3 udir3(inout rng_state_t rng_state) {
	vec2 z = rand2(rng_state);
	vec2 t = vec2(TAU * z.x, acos(2.0 * z.y - 1.0));
	vec2 c = cos(t), s = sin(t);
	return vec3(c.x * s.y, s.x * s.y, c.y);
}

// Box-Muller Transform
// https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
vec2 nrand2(float sigma, vec2 mean, inout rng_state_t rng_state) {
	vec2 z = rand2(rng_state);
	return mean + sigma * sqrt(-2.0 * log(z.x   )) * vec2(cos(TAU * z.y), sin(TAU * z.y));
}

vec3 nrand3(float sigma, vec3 mean, inout rng_state_t rng_state) {
	vec4 z = rand4(rng_state);
	return mean + sigma * sqrt(-2.0 * log(z.xxy )) * vec3(cos(TAU * z.z), sin(TAU * z.z), cos(TAU * z.w));
}

vec4 nrand4(float sigma, vec4 mean, inout rng_state_t rng_state) {
	vec4 z = rand4(rng_state);
	return mean + sigma * sqrt(-2.0 * log(z.xxyy)) * vec4(cos(TAU * z.z), sin(TAU * z.z), cos(TAU * z.w), sin(TAU * z.w));
}

vec2 pixel_filter(vec2 coord, float width, rng_state_t rng_state) {
	return coord + width * (rand2(rng_state) - 0.5); // box filter, uniform over entire pixel
	//return nrand2(width, coord, rng_state); // normal filter, normal distribution
}
*/
