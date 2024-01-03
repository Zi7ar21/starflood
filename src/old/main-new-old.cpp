#include <cstdio>
#include <cstdlib>

// OpenMP
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

typedef float real;

#define inv_tau 0.1591549430918953357688837633725143620344596457404564487476673441

double gaussian2d(double x0, double x1, double sigma) {
	sigma = sigma > 0.0 ? sigma : 1.0;
	double sigma0 = 1.0 / sigma;
	double sigma1 = 1.0 / (2.0*sigma*sigma);
	return sigma0*inv_tau*exp(-sigma1*((x0*x0)+(x1*x1)));
}

int main(int argc, char** argv) {
	if(argc > 1) {
		printf("%s: Invalid Parameter '%s'\n", argv[0], argv[1]);

		return EXIT_FAILURE;
	}

	int w, h;

	double t0, t1;

	t0 = omp_get_wtime();
	t1 = omp_get_wtime();

	size_t image_size = (w * h) * 3 * sizeof(float);

	float* image = (float*)malloc(image_size);

	stbi_flip_vertically_on_write(1);

	const int num_steps = 10;

	const int n = 10; // number of bodies in the simulation

	size_t x_size = n * 3 * sizeof(real);
	size_t v_size = n * 3 * sizeof(real);
	size_t a_size = n * 3 * sizeof(real);
	size_t t_size = n * 3 * sizeof(real);

	real* x = (real*)malloc(x_size); // position buffer
	real* v = (real*)malloc(v_size); // velocity buffer
	real* a = (real*)malloc(a_size); // acceleration buffer
	real* t = (real*)malloc(t_size); // temporary buffer

	for(int i = 0; i < n * 3; i++) x[i] = (real)0;
	for(int i = 0; i < n * 3; i++) v[i] = (real)0;
	for(int i = 0; i < n * 3; i++) a[i] = (real)0;
	for(int i = 0; i < n * 3; i++) t[i] = (real)0;

	const real dt = 0.1;

	for(int step_num; step_num < num_steps; step_num++) {
		// render and write
		{
			// clear the image buffer
			//#pragma omp simd
			for(int i = 0; i < (w * h) * 3; i++) image[i] = (float)0;

			char file_name[FILENAME_MAX];

			snprintf(file_name, FILENAME_MAX, "out/step_%03d.hdr", step_num);

			t0 = omp_get_wtime();

			stbi_write_hdr(file_name, w, h, 3, image);

			t1 = omp_get_wtime();

			//fprintf(logfile, "%.6f,", 1000.0*(t1-t0));
		}

		// kick
		{
			for(int i = 0; i < (n * 3); i++) v[i] += (real)0.5 * dt * a[i];
		}

		// drift
		{
			for(int i = 0; i < (n * 3); i++) x[i] += dt * v[i];
		}

		// update acceleration
		{
			for(int i = 0; i < n * 3; i++) a[i] = (real)0;
		}

		// kick
		{
			for(int i = 0; i < (n * 3); i++) v[i] += (real)0.5 * dt * a[i];
		}
	}

	free(image);

	return EXIT_SUCCESS;
}
