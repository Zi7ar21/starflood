// SPDX-License-Identifier: GPL-3.0-or-later
/* Copyright (c) 2023-2026 Jacob Bingham
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <stdio.h>

#include "config.h"

typedef struct {
	FILE* file;
} log_t;

// Initialize a log
int log_init(log_t* restrict log, const char* restrict filename);

// Free a log
int log_free(log_t* restrict log);

// Sync a log (fflush)
int log_sync(log_t* restrict log);
