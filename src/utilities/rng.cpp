#include <rng.hpp>

/*
// https://nullprogram.com/blog/2018/07/31/
uint32_t triple32(uint32_t x) {
	// exact bias: 0.020888578919738908
	x ^= x >> (uint32_t)17u;
	x *= (uint32_t)0xED5AD4BBu;
	x ^= x >> (uint32_t)11u;
	x *= (uint32_t)0xAC4C1B51u;
	x ^= x >> (uint32_t)15u;
	x *= (uint32_t)0x31848BABu;
	x ^= x >> (uint32_t)14u;
	return x;
}
*/

// http://www.jcgt.org/published/0009/03/02/
rng_state_t pcg4d(rng_state_t v) {
	v.z0 = (uint32_t)0x0019660Du * v.z0 + (uint32_t)0x3C6EF35Fu;
	v.z1 = (uint32_t)0x0019660Du * v.z1 + (uint32_t)0x3C6EF35Fu;
	v.z2 = (uint32_t)0x0019660Du * v.z2 + (uint32_t)0x3C6EF35Fu;
	v.z3 = (uint32_t)0x0019660Du * v.z3 + (uint32_t)0x3C6EF35Fu;

	v.z0 += v.z1*v.z3;
	v.z1 += v.z2*v.z0;
	v.z2 += v.z0*v.z1;
	v.z3 += v.z1*v.z2;

	v.z0 ^= v.z0 >> (uint32_t)16u;
	v.z1 ^= v.z1 >> (uint32_t)16u;
	v.z2 ^= v.z2 >> (uint32_t)16u;
	v.z3 ^= v.z3 >> (uint32_t)16u;

	v.z0 += v.z1*v.z3;
	v.z1 += v.z2*v.z0;
	v.z2 += v.z0*v.z1;
	v.z3 += v.z1*v.z2;

	return v;
}

void update_rng(rng_state_t* rng_state) {
	rng_state_t _rng_state = *rng_state;

	/*
	_rng_state.z0 = triple32(_rng_state.z0);
	_rng_state.z1 = triple32(_rng_state.z1);
	_rng_state.z2 = triple32(_rng_state.z2);
	_rng_state.z3 = triple32(_rng_state.z3);
	*/

	_rng_state = pcg4d(_rng_state);

	*rng_state = _rng_state;
}
