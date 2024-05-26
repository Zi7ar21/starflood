#pragma once

#include <stdint.h>

typedef struct {
	uint32_t s0;
	uint32_t s1;
	uint32_t s2;
	uint32_t s3;
} rng_state_t;

// https://www.jcgt.org/published/0009/03/02/
void update_rng(rng_state_t* rng_state) {
	v = v * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
	v ^= v >> (uint32_t)16u;
	v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
	return v;
}

float urand1(rng_state_t* rng_state) {
}

/*



float urand1(inout rng_state_t rng_state) {
	rng_state.s0 = pcg4d(rng_state.s0);

	return 2.32830629776081821092e-10*float(rng_state.s0.x);
}

vec2 urand2(inout rng_state_t rng_state) {
	rng_state.s0 = pcg4d(rng_state.s0);

	return 2.32830629776081821092e-10*vec2(rng_state.s0.xy);
}

vec3 urand3(inout rng_state_t rng_state) {
	rng_state.s0 = pcg4d(rng_state.s0);

	return 2.32830629776081821092e-10*vec3(rng_state.s0.xyz);
}

vec4 urand4(inout rng_state_t rng_state) {
	rng_state.s0 = pcg4d(rng_state.s0);

	return 2.32830629776081821092e-10*vec4(rng_state.s0);
}

float urand1(rng_state_t* rng_state);

// If defined, PCG will be used for Random Number Generation
// http://www.pcg-random.org/
//#define USE_PCG_RANDOM

// ##### Constants #####


https://float.exposed/
http://www.mimirgames.com/articles/programming/digits-of-pi-needed-for-floating-point-numbers/



It's convenient to make constant floating-point variables with the best possible precision so you don't need to type them everywhere,
and precomputing their inverses allows you to change division (which is slow) into multiplication e.g. x / PI = INV_PI * x


//const float     PI = glm::uintBitsToFloat(0x40490FDBU); // Pi
//const float INV_PI = glm::uintBitsToFloat(0x3EA2F983U); // 1 / Pi
//const float     TAU = glm::uintBitsToFloat(0x40C90FDBU); // 2 * Pi = "Tau"
//const float INV_TAU = glm::uintBitsToFloat(0x3E22F983U); // 1 / (2 * Pi)

const float     PI = 3.1415926535897932384626433832795028841971693993751058209749445923; // Pi
const float INV_PI = 0.3183098861837906715377675267450287240689192914809128974953346881; // 1 / Pi

const float     TAU = 6.2831853071795864769252867665590057683943387987502116419498891846; // 2 * Pi = "Tau"
const float INV_TAU = 0.1591549430918953357688837633725143620344596457404564487476673441; // 1 / (2 * Pi)

// 2^32 (but slightly different, so 32-bit RNG doesn't ever reach 0.0 or 1.0)
// I know I said "C++11 lacks a real way to do this" earlier but there's no real way to get the closest
// but still *slightly* different float besides using hexadecimal
const float     RNG32_MAX = uintBitsToFloat(0x4F7FFFFFU); // Actual: 0x4F800000U
const float INV_RNG32_MAX = uintBitsToFloat(0x2F7FFFFFU); // Actual: 0x2F800000U

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
struct rng_state_t {
	uint seed;
};

void init_rng(inout rng_state_t rng_state, vec2 pixel_coord, vec2 res, int frame_number) {
	pixel_coord = floor(pixel_coord);

	rng_state.seed = uint(frame_number) * uint(res.x * res.y) + uint(pixel_coord.x + res.x * pixel_coord.y) + 1u;
}

// Triple32 Hash: https://nullprogram.com/blog/2018/07/31/
// this hash in particular is remarkable because it's a statistically perfect 32-bit integer hash (of this kind)
uint triple32(uint x) {
	x ^= x >> 17u;
	x *= 0xED5AD4BBu;
	x ^= x >> 11u;
	x *= 0xAC4C1B51u;
	x ^= x >> 15u;
	x *= 0x31848BABu;
	x ^= x >> 14u;
	return x;
}

float rand(inout rng_state_t rng_state) {
	rng_state.seed = triple32(rng_state.seed);

	return INV_RNG32_MAX * float(rng_state.seed);
}

vec2 rand2(inout rng_state_t rng_state) {
	vec2 v;
	v.x = rand(rng_state);
	v.y = rand(rng_state);
	return v;
}

vec3 rand3(inout rng_state_t rng_state) {
	vec3 v;
	v.x = rand(rng_state);
	v.y = rand(rng_state);
	v.z = rand(rng_state);
	return v;
}

vec4 rand4(inout rng_state_t rng_state) {
	vec4 v;
	v.x = rand(rng_state);
	v.y = rand(rng_state);
	v.z = rand(rng_state);
	v.w = rand(rng_state);
	return v;
}
#endif

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