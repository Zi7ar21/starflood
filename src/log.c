#include "log.h"

#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "common.h"

int log_init(log_t* restrict log, const char* restrict filename) {
	if(NULL == (void*)filename) {
		return STARFLOOD_FAILURE;
	}

	log->file = fopen(filename, "w");

	if(NULL == (void*)log->file) {
		fprintf(stderr, "%s error: fopen(\"%s\", \"%s\") ", "log_init()", filename, "w");

		perror("failed");

		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}

int log_free(log_t* restrict log) {
	if( 0 != fclose(log->file) ) {
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}

int log_sync(log_t* restrict log) {
	if( 0 != fflush(log->file) ) {
		return STARFLOOD_FAILURE;
	}

	return STARFLOOD_SUCCESS;
}
