// SPDX-License-Identifier: GPL-3.0-or-later
/* Copyright (c) 2023-2026 Jacob Bingham
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <stdint.h>

#include "pcg.h"

// The closest (but still less than) IEEE-754 binary32 value to 1 / exp2(32)
// Multiplying a random, uniformly distributed unsigned 32-bit integer by
// this value should produce a random real number in the range [0.0, 1.0).
// See: https://float.exposed/0x2F7FFFFF
#define INV_PCG32_MAX 2.3283062977608182109179324470460414886474609375e-10

// 3-round multiply-xorshift 1-dimensional hash function
// Bias: ~0.020888578919738908 (theoretical limit for a hash of this kind)
// Note: Chris Wellons has since updated his article with links explaining
// caveats to this kind of hashing.
// See: https://nullprogram.com/blog/2018/07/31/
static inline void triple32(uint32_t* restrict s) {
	uint32_t x = *s;
	x ^= x >> (uint32_t)17u;
	x *= (uint32_t)0xED5AD4BBu;
	x ^= x >> (uint32_t)11u;
	x *= (uint32_t)0xAC4C1B51u;
	x ^= x >> (uint32_t)15u;
	x *= (uint32_t)0x31848BABu;
	x ^= x >> (uint32_t)14u;
	*s = x;
}
