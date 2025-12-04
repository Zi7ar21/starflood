#pragma once

#include <stddef.h>

#include "types.h"

// Memory alignment
// Default: 4 MiB (4 * 1024 * 1024 bytes)
#define STARFLOOD_ALIGNMENT (size_t)4194304u

// Gravitational constant (physical units are not yet implemented)
#define G (real)1.0

// alternative (m^3 * kg^-1 * s^2)
// https://physics.nist.gov/cgi-bin/cuu/Value?bg
//#define G 6.6743e-11
