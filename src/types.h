#pragma once

#include <stdint.h>

#include "config.h"

typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef  float f32;
typedef double f64;

// Floating-point type used for simulation
#ifdef STARFLOOD_DOUBLE_PRECISION
typedef f64 real;
#else
typedef f32 real;
#endif
