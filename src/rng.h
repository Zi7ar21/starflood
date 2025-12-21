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
inline void pcg4d(uint32_t* restrict s) {
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

/*
#ifdef _OPENMP
	#ifdef ENABLE_OFFLOADING
	#pragma omp end declare target
	#endif
#endif
*/
