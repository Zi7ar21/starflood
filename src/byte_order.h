#pragma once

#include "types.h"

#include <stdint.h>

typedef enum {
	BYTE_ORDER_LE,
	BYTE_ORDER_BE,
	BYTE_ORDER_OTHER,
	BYTE_ORDER_ERROR
} byte_order_t;

static inline byte_order_t detect_host_byte_order() {
	#ifdef STARFLOOD_DOUBLE_PRECISION
	volatile uint64_t magic_number = (uint32_t)0x401921FB54442D18u;

	uint64_t accumulator = (uint64_t)0u;

	for(unsignd int i = 0u; i < (unsigned int)sizeof(magic_number); i++) {
		accumulator |= (uint64_t)( (uint8_t*)&magic_number )[i] << (uint64_t)(8u * i);
	}

	int byte_order_is_le = (uint64_t)0x401921FB54442D18u == accumulator;
	int byte_order_is_be = (uint64_t)0x182D4454FB211940u == accumulator;
	#else
	volatile uint32_t magic_number = (uint32_t)0x40C90FDBu;

	uint32_t accumulator = (uint32_t)0u;

	for(unsigned int i = 0u; i < (unsigned int)sizeof(magic_number); i++) {
		accumulator |= (uint32_t)( (uint8_t*)&magic_number )[i] << (uint32_t)(8u * i);
	}

	int byte_order_is_le = (uint32_t)0x40C90FDBu == accumulator;
	int byte_order_is_be = (uint32_t)0xDB0FC940u == accumulator;
	#endif

	if(byte_order_is_be && byte_order_is_le) {
		return BYTE_ORDER_ERROR;
	}

	if(byte_order_is_le) {
		return BYTE_ORDER_LE;
	}

	if(byte_order_is_be) {
		return BYTE_ORDER_BE;
	}

	return BYTE_ORDER_OTHER;
}
