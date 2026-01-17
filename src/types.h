#pragma once

#include "config.h"

#include <stdint.h>

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

// Unsigned integer types with same size as real type (for alignment purposes)
#ifdef STARFLOOD_DOUBLE_PRECISION
typedef u64 uint;
#else
typedef u32 uint;
#endif

// String of the floating-point type name (used for writing .ply file header)
#ifdef STARFLOOD_DOUBLE_PRECISION
#define REAL_TYPE_STRING "double"
#else
#define REAL_TYPE_STRING "float"
#endif

// floating-point math function aliases
#ifdef STARFLOOD_DOUBLE_PRECISION
#define real_cos(__x) (real)cos(__x)
#define real_sin(__x) (real)sin(__x)
#define real_tan(__x) (real)sin(__x)
#define real_exp(__x) (real)exp(__x)
#define real_exp2(__x) (real)exp2(__x)
#define real_log(__x) (real)log(__x)
#define real_log2(__x) (real)log2(__x)
#define real_cbrt(__x) (real)cbrt(__x)
#define real_sqrt(__x) (real)sqrt(__x)
#define real_pow(__x,__y) (real)pow((__x),(__y))
#else
#define real_cos(__x) (real)cosf(__x)
#define real_sin(__x) (real)sinf(__x)
#define real_tan(__x) (real)sinf(__x)
#define real_exp(__x) (real)expf(__x)
#define real_exp2(__x) (real)exp2f(__x)
#define real_log(__x) (real)logf(__x)
#define real_log2(__x) (real)log2f(__x)
#define real_cbrt(__x) (real)cbrtf(__x)
#define real_sqrt(__x) (real)sqrtf(__x)
#define real_pow(__x,__y) (real)powf((__x),(__y))
#endif
