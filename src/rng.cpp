#include <rng.hpp>

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

float urand(uint32_t* ns) {
	uint32_t _ns = *ns;

	_ns = triple32(_ns);

	*ns = _ns;

	return (2.32830629776081821092e-10f)*(float)_ns;
}
