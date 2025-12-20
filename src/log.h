#pragma once

#include <stdio.h>

#include "config.h"

typedef struct {
	FILE* file;
} log_t;

int log_init(log_t* restrict log, const char* restrict filename);

int log_free(log_t* restrict log);

int log_sync(log_t* restrict log);
