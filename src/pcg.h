// SPDX-License-Identifier: MIT
/*  Copyright (c) 2020-2026 Mark Jarzynski, Marc Olano, and Jacob Bingham
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#pragma once

#include <stdint.h>

// Permuted congruential generator 4-dimensional hash function
// See: https://www.jcgt.org/published/0009/03/02/
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
