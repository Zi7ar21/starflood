// SPDX-License-Identifier: GPL-3.0-or-later
/* Copyright (c) 2023-2026 Jacob Bingham
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "simulation.h"

#ifdef ENABLE_SPH
int solve_sph(const sim_t* restrict sim_ptr);
#endif
