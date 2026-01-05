#pragma once

#include "config.h"

#ifndef ENABLE_TIMING
	#define TIMING_INIT()
	#define TIMING_START()
	#define TIMING_STOP()
	#define TIMING_PRINT(function_name,label_name) 
	#define LOG_TIMING(log)
#else
	#ifndef TIMING_USE_OMP_GET_WTIME
		#if _POSIX_C_SOURCE >= 199309L
			#include <stdint.h>
			#include <time.h>

			#ifndef STARFLOOD_POSIX_CLOCKID
				#define STARFLOOD_POSIX_CLOCKID CLOCK_MONOTONIC
			#endif

			#define TIMING_INIT() struct timespec __t0, __t1;

			#define TIMING_START() \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t1); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t0); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t1); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t0); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t1); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t0); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t1); \
			clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t0)

			#define TIMING_STOP() clock_gettime(STARFLOOD_POSIX_CLOCKID, &__t1)

			#define TIMING_PRINT(function_name,label_name) printf("%s: %.03f usec %s\n",function_name,1.0e-3*(double)((intmax_t)1000000000l*(intmax_t)(__t1.tv_sec-__t0.tv_sec)+((intmax_t)__t1.tv_nsec-(intmax_t)__t0.tv_nsec)),label_name); fflush(stdout)

			#define LOG_TIMING(log) fprintf((log).file,",%.03f",1.0e-3*(double)((intmax_t)1000000000l*(intmax_t)(__t1.tv_sec-__t0.tv_sec)+((intmax_t)__t1.tv_nsec-(intmax_t)__t0.tv_nsec)))
		#else
			#error "When not using #define TIMING_USE_OMP_GET_WTIME, POSIX function clock_gettime() is required, which needs _POSIX_C_SOURCE >= 199309L to be defined before including headers from the standard library."
		#endif
	#else 
		#ifndef _OPENMP
			#error "When using #define TIMING_USE_OMP_GET_WTIME, OpenMP function omp_get_wtime() is required, which needs _OPENMP to be defined. (Are you sure the compiler has OpenMP enabled?)"
		#endif

		#define TIMING_INIT() \
		volatile double __t0 = 0.0; \
		volatile double __t1 = 0.0

		#define TIMING_START() \
		__t1 = omp_get_wtime(); \
		__t0 = omp_get_wtime(); \
		__t1 = omp_get_wtime(); \
		__t0 = omp_get_wtime(); \
		__t1 = omp_get_wtime(); \
		__t0 = omp_get_wtime(); \
		__t1 = omp_get_wtime(); \
		__t0 = omp_get_wtime()

		#define TIMING_STOP() __t1 = omp_get_wtime()

		#define TIMING_PRINT(function_name,label_name) printf("%s: %.06f usec %s\n",function_name,1.0e6*(__t1-__t0),label_name)

		#define LOG_TIMING(log) fprintf((log).file,"%.03f",1.0e6*(__t1-__t0))
	#endif
#endif
