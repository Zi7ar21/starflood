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
