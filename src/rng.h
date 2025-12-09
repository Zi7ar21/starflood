#pragma once

#include <stdint.h>

// https://float.exposed/0x2F7FFFFF
#define INV_PCG32_MAX 2.3283062977608182109179324470460414886474609375e-10

// https://www.jcgt.org/published/0009/03/02/
void pcg4d(uint32_t* restrict s);
