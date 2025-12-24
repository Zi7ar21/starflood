#pragma once

#include <stdint.h>

// https://float.exposed/0x2F7FFFFF
#define INV_PCG32_MAX 2.3283062977608182109179324470460414886474609375e-10

/*
#ifdef _OPENMP
	#ifdef ENABLE_OFFLOADING
	#pragma omp declare target
	#endif
#endif
*/

// https://www.jcgt.org/published/0009/03/02/
static inline void pcg4d(uint32_t* restrict s) {
	uint32_t v[4] = {s[0], s[1], s[2], s[3]};

	v[0] = v[0] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v[1] = v[1] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v[2] = v[2] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v[3] = v[3] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;

	v[0] += v[1] * v[3];
	v[1] += v[2] * v[0];
	v[2] += v[0] * v[1];
	v[3] += v[1] * v[2];

	v[0] = v[0] ^ (v[0] >> (uint32_t)16u);
	v[1] = v[1] ^ (v[1] >> (uint32_t)16u);
	v[2] = v[2] ^ (v[2] >> (uint32_t)16u);
	v[3] = v[3] ^ (v[3] >> (uint32_t)16u);

	v[0] += v[1] * v[3];
	v[1] += v[2] * v[0];
	v[2] += v[0] * v[1];
	v[3] += v[1] * v[2];

	s[0] = v[0];
	s[1] = v[1];
	s[2] = v[2];
	s[3] = v[3];
}

// https://nullprogram.com/blog/2018/07/31/
static inline void triple32(uint32_t* restrict s) {
	uint32_t x = *s;

	x ^= x >> (uint32_t)17u;
	x *= (uint32_t)0xED5AD4BBu;
	x ^= x >> (uint32_t)11u;
	x *= (uint32_t)0xAC4C1B51u;
	x ^= x >> (uint32_t)15u;
	x *= (uint32_t)0x31848BABu;
	x ^= x >> (uint32_t)14u;

	*s = x;
}

/*
#ifdef _OPENMP
	#ifdef ENABLE_OFFLOADING
	#pragma omp end declare target
	#endif
#endif
*/
