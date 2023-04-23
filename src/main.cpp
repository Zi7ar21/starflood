#include <iostream>
#include <cstdint>
//#include <cstdlib>
#include <cmath>

//#include <Eigen/Core>

#include "common.hpp"
#include "config.h"
#include "constants.h"
#include "image.hpp"

/*
// exact bias: 0.020888578919738908
uint32_t triple32(uint32_t x) {
	x ^= x >> (uint32_t)17U;
	x *= (uint32_t)0xED5AD4BBU;
	x ^= x >> (uint32_t)11U;
	x *= (uint32_t)0xAC4C1B51U;
	x ^= x >> (uint32_t)15U;
	x *= (uint32_t)0x31848BABU;
	x ^= x >> (uint32_t)14U;
	return x;
}
*/


/*
def W( x, y, z, h ):
	"""
	Gausssian Smoothing kernel (3D)
	x     is a vector/matrix of x positions
	y     is a vector/matrix of y positions
	z     is a vector/matrix of z positions
	h     is the smoothing length
	w     is the evaluated smoothing function
	"""
	
	r = np.sqrt(x**2 + y**2 + z**2)
	
	w = (1.0 / (h*np.sqrt(np.pi)))**3 * np.exp( -r**2 / h**2)
	
	return w
*/

//float W(float x, float y, float z) {
//}

void UpdateAcceleration(unsigned int N, float* m, float* r, float* v, float* a) {
	// Compute Forces (O(N^2), very slow...)
	{
		for(unsigned int i = 0; i < N; i++) {
			float m1 = m[i];

			float xi[2] = {{r[2 * i + 0]}, {r[2 * i + 1]}};

			//#pragma omp parallel for
			for(unsigned int j = 0; j < N; j++) {
				if(i == j) continue;

				float m2 = m[j];

				float d[2] = {
					{r[2 * j + 0] - xi[0]},
					{r[2 * j + 1] - xi[1]},
				};

				float r2 = (d[0] * d[0]) + (d[1] * d[1]);
				float r1  = sqrtf(r2);

				a[2 * i + 0] += d[0] * G_SI * m1 * m2 * ((float)1 / (sqrtf(r1) * r2));
				a[2 * i + 1] += d[1] * G_SI * m1 * m2 * ((float)1 / (sqrtf(r1) * r2));
			}

			a[2 * i + 0] += -0.001 * xi[0];
			a[2 * i + 1] += -0.001 * xi[1];
		}
	}
}

int main(int argc, char** argv) {
	if(argc > 1) {
		std::cerr << "Invalid argument: \"" << argv[1] << "\"" << std::endl;

		return EXIT_FAILURE;
	}

	unsigned int N = 2048; // number of particles

	double dt = 0.06666666666666666666666666666667;
	//double dt = 0.66666666666666666666666666666667;
	//double dt = 1.0;

	unsigned int ImageSizeX = 512;
	unsigned int ImageSizeY = 512;

	std::cout
	<< "+---------------+\n"
	<< "|   Starflood   |\n"
	<< "+---------------+\n"
	<< "\n"
	<< "Parameters:\n"
	<< "{\n"
	<< "  N: " << N << "\n"
	<< "  dt: " << dt << "\n"
	<< "  ImageSizeX: " << ImageSizeX << "\n"
	<< "  ImageSizeY: " << ImageSizeY << "\n"
	<< "}\n"
	<< std::endl;

	//return EXIT_SUCCESS;

	Buffer2D<RGBA32F> image(ImageSizeX, ImageSizeY);

	//simulation.Initialize();

	// row-major
	float* m = (float*)malloc(sizeof(float) * (size_t)N * (size_t)1); // Particle masses
	float* r = (float*)malloc(sizeof(float) * (size_t)N * (size_t)2); // Particle positions
	float* v = (float*)malloc(sizeof(float) * (size_t)N * (size_t)2); // Particle velocities
	float* a = (float*)malloc(sizeof(float) * (size_t)N * (size_t)2); // Particle accelerations

	// initialize particle mass
	//for(unsigned int i = 0; i < N; i++) m[i] = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
	for(unsigned int i = 0; i < N; i++) m[i] = 100000.0/N;

	std::srand(1);

	// initialize particle position
	for(unsigned int i = 0; i < (2 * N); i++) r[i] = 1.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);

	std::srand(2);

	// initialize particle velocity
	for(unsigned int i = 0; i < (2 * N); i++) v[i] = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
	//for(unsigned int i = 0; i < (2 * N); i++) v[i] = 0.0;
	//for(unsigned int i = 0; i < (2 * N); i++) {
	//	v[2 * i + 0] =  0.01 * r[2 * i + 1] / sqrtf(r[2 * i + 0]*r[2 * i + 0]+r[2 * i + 1]*r[2 * i + 1]);
	//	v[2 * i + 1] = -0.01 * r[2 * i + 0] / sqrtf(r[2 * i + 0]*r[2 * i + 0]+r[2 * i + 1]*r[2 * i + 1]);
	//}

	UpdateAcceleration(N, m, r, v, a); // update acceleration

	double t = 0.0;

	for(int n = 1; n < 300; n++) {
		//#pragma omp simd
		for(unsigned int i = 0; i < (ImageSizeX * ImageSizeY); i++) image.data[i] = RGBA32F(0.0, 0.0, 0.0, 1.0); // Clear image buffer

		// Render the image
		for(unsigned int i = 0; i < N; i++) {
			// vec2 uv = (fragCoord - 0.5 * iResolution.xy) / min(iResolution.x, iResolution.y);
			// vec2 fragCoord = (min(iResolution.x, iResolution.y) * uv) + (0.5 * iResolution.xy);

			//Eigen::Vector3f p = r.row(i);

			//double _t = 0.0125 * t;

			//Eigen::Vector3f _p(
			//(Real)(p(0) * cos(2.0 * acos(-1.0) * _t)) - (p(2) * sin(2.0 * acos(-1.0) * _t)),
			//(Real)(p(0) * sin(2.0 * acos(-1.0) * _t)) + (p(2) * sin(2.0 * acos(-1.0) * _t)),
			//p(1)
			//);

			//p = _p;

			float px = r[2 * i + 0];
			float py = r[2 * i + 1];

			float uv_x = (fminf(ImageSizeX, ImageSizeY) * (3e-1 * px)) + (0.5 * ImageSizeX);
			float uv_y = (fminf(ImageSizeX, ImageSizeY) * (3e-1 * py)) + (0.5 * ImageSizeY);

			int coord_x = (int)uv_x;
			int coord_y = (int)uv_y;

			if(0 <= coord_x && coord_x < ImageSizeX && 0 <= coord_y && coord_y < ImageSizeY) {
				image.data[coord_x + ImageSizeX * coord_y].r += 1.0;
				image.data[coord_x + ImageSizeX * coord_y].g += 1.0;
				image.data[coord_x + ImageSizeX * coord_y].b += 1.0;
			}
		}

		// Write the image buffer
		{
			char filename[80];

			sprintf(filename, "out/IMG_%04d.tga", n);
			//sprintf(filename, "out/IMG_%04d.hdr", n);

			WriteImageBuffer(image, filename);
			//WriteImageBufferHDR(image, filename);
		}

		//UpdateSimulation(simulation);

		//#pragma omp simd
		for(unsigned int i = 0; i < (2 * N); i++) v[i] += 0.5 * dt * a[i]; // 1/2 kick
		//v = v + (0.5 * dt * a); // 1/2 kick

		for(unsigned int i = 0; i < (2 * N); i++) v[i] *= 0.99; // drift

		//#pragma omp simd
		for(unsigned int i = 0; i < (2 * N); i++) r[i] += dt * v[i]; // drift
		//r = r + (dt * v); // drift

		UpdateAcceleration(N, m, r, v, a); // update acceleration

		//#pragma omp simd
		for(unsigned int i = 0; i < (2 * N); i++) v[i] += 0.5 * dt * a[i]; // 1/2 kick
		//v = v + (0.5 * dt * a); // 1/2 kick

		for(unsigned int i = 0; i < (2 * N); i++) v[i] *= 0.99; // drift

		t += dt;
	}

	//free(m);
	//free(x);
	//free(v);
	//free(a);

	//simulation.Deallocate();

	std::cout << "Finished!" << std::endl;

	return EXIT_SUCCESS;
}
