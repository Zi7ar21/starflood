#pragma once

#include <stdint.h>

// https://float.exposed/0x2F7FFFFF
#define INV_PCG32_MAX 2.32830629776081821092e-10

// https://www.jcgt.org/published/0009/03/02/
void pcg4d(uint32_t* restrict s);
