#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <math.h>

// OpenMP
#include <omp.h>

#include <fftw3.h>

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

	int w = 255, h = 255;

	double t0, t1;

	t0 = omp_get_wtime();
	t1 = omp_get_wtime();

	size_t image_size = (w * h) * 3 * sizeof(float);

	float* image = (float*)malloc(image_size);

	stbi_flip_vertically_on_write(1);

	size_t mass_i_size   = (w * h) * sizeof(double);
	size_t mass_o_size   = (w * h) * sizeof(fftw_complex);
	size_t kern_x_i_size = (w * h) * sizeof(double);
	size_t kern_x_o_size = (w * h) * sizeof(fftw_complex);
	size_t kern_y_i_size = (w * h) * sizeof(double);
	size_t kern_y_o_size = (w * h) * sizeof(fftw_complex);
	size_t solv_x_i_size = (w * h) * sizeof(fftw_complex);
	size_t solv_x_o_size = (w * h) * sizeof(fftw_complex);
	size_t solv_y_i_size = (w * h) * sizeof(fftw_complex);
	size_t solv_y_o_size = (w * h) * sizeof(fftw_complex);

	fftw_complex* mass_i   = (fftw_complex*)fftw_malloc(mass_i_size);
	fftw_complex* mass_o   = (fftw_complex*)fftw_malloc(mass_o_size);
	fftw_complex* kern_x_i = (fftw_complex*)fftw_malloc(kern_x_i_size);
	fftw_complex* kern_x_o = (fftw_complex*)fftw_malloc(kern_x_o_size);
	fftw_complex* kern_y_i = (fftw_complex*)fftw_malloc(kern_y_i_size);
	fftw_complex* kern_y_o = (fftw_complex*)fftw_malloc(kern_y_o_size);
	fftw_complex* solv_x_i = (fftw_complex*)fftw_malloc(solv_x_i_size);
	fftw_complex* solv_x_o = (fftw_complex*)fftw_malloc(solv_x_o_size);
	fftw_complex* solv_y_i = (fftw_complex*)fftw_malloc(solv_y_i_size);
	fftw_complex* solv_y_o = (fftw_complex*)fftw_malloc(solv_y_o_size);

	for(int i = 0; i < h; i++) {
		for(int j = 0; j < w; j++) {
			int x = j-93;
			int y = i-163;

			mass_i[i*w+j] = (double)0;
			if((x == 0) && (y == 0)) mass_i[i*w+j] = (double)1;

			x = j-127;
			y = i-127;

			kern_x_i[i*w+j] = -((double)x/(double)sqrt(((double)x*(double)x)+((double)y*(double)y)))/(((double)x*(double)x)+((double)y*(double)y));
			kern_y_i[i*w+j] = -((double)y/(double)sqrt(((double)x*(double)x)+((double)y*(double)y)))/(((double)x*(double)x)+((double)y*(double)y));
			if((x == 0) && (y == 0)) kern_x_i[i*w+j] = (double)0;
			if((x == 0) && (y == 0)) kern_y_i[i*w+j] = (double)0;
		}
	}

	fftw_plan  massfft = fftw_plan_dft_r2c_2d(w, h,   mass_i,   mass_o, FFTW_ESTIMATE); fftw_execute( massfft);
	fftw_plan kernxfft = fftw_plan_dft_r2c_2d(w, h, kern_x_i, kern_x_o, FFTW_ESTIMATE); fftw_execute(kernxfft);
	fftw_plan kernyfft = fftw_plan_dft_r2c_2d(w, h, kern_y_i, kern_y_o, FFTW_ESTIMATE); fftw_execute(kernyfft);

	fftw_destroy_plan( massfft);
	fftw_destroy_plan(kernxfft);
	fftw_destroy_plan(kernyfft);

	for(int i = 0; i < (w * h); i++) {
			solv_x_i[i][0] = mass_o[i][0]*kern_x_o[i][0]-mass_o[i][1]*kern_x_o[i][1];
			solv_x_i[i][1] = mass_o[i][0]*kern_x_o[i][1]+mass_o[i][1]*kern_x_o[i][0];

			solv_y_i[i][0] = mass_o[i][0]*kern_y_o[i][0]-mass_o[i][1]*kern_y_o[i][1];
			solv_y_i[i][1] = mass_o[i][0]*kern_y_o[i][1]+mass_o[i][1]*kern_y_o[i][0];
	}

	fftw_plan solvxfft = fftw_plan_dft_2d(w, h, solv_x_i, solv_x_o, FFTW_BACKWARD, FFTW_ESTIMATE); fftw_execute(solvxfft);
	fftw_plan solvyfft = fftw_plan_dft_2d(w, h, solv_y_i, solv_y_o, FFTW_BACKWARD, FFTW_ESTIMATE); fftw_execute(solvyfft);

	fftw_destroy_plan(solvxfft);
	fftw_destroy_plan(solvyfft);

	// render and write
	{
		// clear the image buffer
		//#pragma omp simd
		for(int i = 0; i < (w * h); i++) {
			image[3*i+0] = 0.1*sqrt(solv_x_o[i][0]*solv_x_o[i][0]+solv_x_o[i][1]*solv_x_o[i][1]);
			image[3*i+1] = 0.1*sqrt(solv_y_o[i][0]*solv_y_o[i][0]+solv_y_o[i][1]*solv_y_o[i][1]);
			image[3*i+2] = (float)0;

			//image[3*i+0] = kern_x_o[i][0];
			//image[3*i+1] = kern_x_o[i][1];
			//image[3*i+2] = 0.1*sqrt(image[3*i+0]*image[3*i+0]+image[3*i+1]*image[3*i+1]);

			//image[3*i+0] = kern_x_i[i];
			//image[3*i+1] = kern_y_i[i];
			//image[3*i+2] = mass_i[i];

			image[3*i+0] = fmaxf(image[3*i+0],(float)0);
			image[3*i+1] = fmaxf(image[3*i+1],(float)0);
			image[3*i+2] = fmaxf(image[3*i+2],(float)0);
		}

		char file_name[FILENAME_MAX];

		snprintf(file_name, FILENAME_MAX, "out/step_%03d.hdr", 0);

		t0 = omp_get_wtime();

		stbi_write_hdr(file_name, w, h, 3, image);

		t1 = omp_get_wtime();

		//fprintf(logfile, "%.9f,", 1000.0*(t1-t0));
	}

	fftw_free(mass_i);
	fftw_free(mass_o);
	fftw_free(kern_x_i);
	fftw_free(kern_x_o);
	fftw_free(kern_y_i);
	fftw_free(kern_y_o);
	fftw_free(solv_x_i);
	fftw_free(solv_x_o);
	fftw_free(solv_y_i);
	fftw_free(solv_y_o);

	/*
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

		//fprintf(logfile, "%.9f,", 1000.0*(t1-t0));
	}
	*/

	free(image);

	return EXIT_SUCCESS;
}
