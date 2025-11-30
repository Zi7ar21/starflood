#ifdef _OPENMP
#pragma message("OpenMP is ENABLED")
#include <omp.h>
#else
#pragma message("OpenMP is DISABLED")
#endif

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

void pcg4d(uint32_t* restrict s) {
	uint32_t v[4] = {s[0], s[1], s[2], s[3]};

	v[0] = v[0] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v[1] = v[1] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v[2] = v[2] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;
	v[3] = v[3] * (uint32_t)0x0019660Du + (uint32_t)0x3C6EF35Fu;

	v[0] += v[1]*v[3];
	v[1] += v[2]*v[0];
	v[2] += v[0]*v[1];
	v[3] += v[1]*v[2];

	v[0] = v[0] ^ (v[0] >> (uint32_t)16u);
	v[1] = v[1] ^ (v[1] >> (uint32_t)16u);
	v[2] = v[2] ^ (v[2] >> (uint32_t)16u);
	v[3] = v[3] ^ (v[3] >> (uint32_t)16u);

	v[0] += v[1]*v[3];
	v[1] += v[2]*v[0];
	v[2] += v[0]*v[1];
	v[3] += v[1]*v[2];

	s[0] = v[0];
	s[1] = v[1];
	s[2] = v[2];
	s[3] = v[3];
}

typedef float real;

int main(int argc, char *argv[]) {
	printf("=== starflood ===\n");

	fflush(stdout);

	#ifdef _OPENMP
	printf("_OPENMP: %d\n", _OPENMP);
	#endif

	unsigned int num_timesteps = 900u;

	if(argc > 1) {
		fprintf(stderr, "%s error: unrecognized option \"%s\"\n", argv[0], argv[1]);

		return EXIT_FAILURE;
	}

	#ifdef _OPENMP
	int max_threads = omp_get_max_threads();

	printf("max_threads: %d\n", max_threads);
	#endif

	//printf("");

	int rendering_enabled = true;

	// Number of bodies to simulate
	size_t N = (size_t)10000u;

	size_t u_size = N * (size_t)1u;
	size_t m_size = N * (size_t)1u;
	size_t r_size = N * (size_t)3u;
	size_t v_size = N * (size_t)3u;
	size_t a_size = N * (size_t)3u;

	size_t u_offset = (size_t)0u;
	size_t m_offset = u_offset + u_size;
	size_t r_offset = m_offset + m_size;
	size_t v_offset = r_offset + r_size;
	size_t a_offset = v_offset + v_size;
	size_t buf_size = a_offset + a_size;

	size_t shared_memory_size = sizeof(real)*buf_size;

	const size_t align_size = (size_t)4096u; // Set to page size or something similar

	void* shared_memory = aligned_alloc(align_size, shared_memory_size);

	if(NULL == shared_memory) {
		perror("error allocating shared_memory");

		return EXIT_FAILURE;
	}

	printf("shared_memory: %p\n", shared_memory);

	memset(shared_memory, 0, shared_memory_size);

	real* buffer = (real*)shared_memory;

	real* u = &(buffer[u_offset]);
	real* m = &(buffer[m_offset]);
	real* r = &(buffer[r_offset]);
	real* v = &(buffer[v_offset]);
	real* a = &(buffer[a_offset]);

	for(size_t i = (size_t)0u; i < N; i++) {
		u[i] = (real)0.0;
	}

	for(size_t i = (size_t)0u; i < N; i++) {
		m[i] = (double)1.0/(double)N;
	}

	for(size_t i = (size_t)0u; i < N; i++) {
		uint32_t s[4] = {(uint32_t)i, (uint32_t)420u, (uint32_t)69u, (uint32_t)1337u};

		pcg4d(s);

		r[3*i+0] = 2.0*((double)s[0]/(double)0xFFFFFFFFu)-1.0;
		r[3*i+1] = 2.0*((double)s[1]/(double)0xFFFFFFFFu)-1.0;
		r[3*i+2] = 2.0*((double)s[2]/(double)0xFFFFFFFFu)-1.0;
	}

	for(size_t i = (size_t)0u; i < (size_t)3u*N; i++) {
		v[i] = (real)0.0;
	}

	for(size_t i = (size_t)0u; i < (size_t)3u*N; i++) {
		a[i] = (real)0.0;
	}

	const real timestep = (real)0.314159;
	const real G = (real)1.0;

	size_t image_w = (size_t)480u;
	size_t image_h = (size_t)480u;

	float* image = (float*)aligned_alloc(align_size, sizeof(float)*(size_t)3*image_w*image_h);

	if(NULL == image) {
		perror("error allocating image");

		free(shared_memory);

		return EXIT_FAILURE;
	}

	printf("image: %p\n", (void*)image);

	for(size_t i = (size_t)0u; i < (size_t)3u*image_w*image_h; i++) {
		image[i] = (float)0.0;
	}

	for(unsigned int n = 0u; n < num_timesteps; n++) {
		if(rendering_enabled) {
			for(size_t y = (size_t)0u; y < image_h; y++) {
				for(size_t x = (size_t)0u; x < image_w; x++) {
					image[3*(image_w*y+x)+0] = (float)0;
					image[3*(image_w*y+x)+1] = (float)0;
					image[3*(image_w*y+x)+2] = (float)0;
				}
			}

			for(unsigned int i = 0u; i < N; i++) {
				int x = ( (float)image_h*((float)0.25*r[3*i+0]) )+((float)0.5*(float)image_w)-(float)0.5;
				int y = ( (float)image_h*((float)0.25*r[3*i+1]) )+((float)0.5*(float)image_h)-(float)0.5;

				if((int)0 <= x && x < (int)image_w && (int)0 <= y && y < (int)image_h) {
					image[3*(image_w*y+x)+0] += (float)0.333;
					image[3*(image_w*y+x)+1] += (float)0.333;
					image[3*(image_w*y+x)+2] += (float)0.333;
				}
			}

			char filename[FILENAME_MAX];

			snprintf(filename, FILENAME_MAX, "./out/%04d.pfm", n);

			FILE* file = fopen(filename, "w");

			if(NULL == file) {
				fprintf(stderr, "fopen(%s, \"w\") failed", filename);

				perror("");

				break;
			}

			fprintf(file, "PF\n%zu %zu\n-1.0\n", image_w, image_h);

			for(size_t y = (size_t)0u; y < image_h; y++) {
				for(size_t x = (size_t)0u; x < image_w; x++) {
					float pixel[3] = {
						image[3*(image_w*y+x)+0],
						image[3*(image_w*y+x)+1],
						image[3*(image_w*y+x)+2]
					};

					fwrite(pixel, sizeof(float), (size_t)3, file);
				}
			}

			fclose(file);
		}

		#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic,100)
		#endif
		for(unsigned int i = 0u; i < N; i++) {
			real U_sum = (real)0;

			real F[3] = {(real)0, (real)0, (real)0};

			for(unsigned int j = 0u; j < N; j++) {
				if(i == j) {
					continue;
				}

				real dr[3] = {r[3*i+0]-r[3*j+0], r[3*i+1]-r[3*j+1], r[3*i+2]-r[3*j+2]};

				real r2 = (dr[0]*dr[0])+(dr[1]*dr[1])+(dr[2]*dr[2]);

				real r1 = (real)sqrtf(r2);

				real rinv = (real)1.0 / r1;

				real 

				U = -G*m[i]*m[j]*rinv;

				U_sum += U;

				F[0] += U*dr[0]/(r2+(real)0.001);
				F[1] += U*dr[1]/(r2+(real)0.001);
				F[2] += U*dr[2]/(r2+(real)0.001);
			}

			u[i] = U_sum;
			a[3*i+0] = F[0];
			a[3*i+1] = F[1];
			a[3*i+2] = F[2];
		}

		double sum = 0.0;

		for(unsigned int i = 0u; i < N; i++) {
			sum += u[i];
		}

		printf("%.015f\n", sum);

		for(unsigned int i = 0u; i < (size_t)3*N; i++) {
			v[i] += timestep*a[i];
		}

		for(unsigned int i = 0u; i < (size_t)3*N; i++) {
			r[i] += timestep*v[i];
		}
	}

	free(image);
	free(shared_memory);

	return EXIT_SUCCESS;
}
