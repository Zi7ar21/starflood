#include <iostream>
#include <cstdint>
//#include <cstdlib>
//#include <cmath>

#include <Eigen/Core>

#include "common.hpp"
#include "config.h"
#include "constants.h"
#include "image.hpp"

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

void UpdateAcceleration(
unsigned int N,
Eigen::Matrix<Real, Eigen::Dynamic, 1>& m,
Eigen::Matrix<Real, Eigen::Dynamic, 3>& r,
Eigen::Matrix<Real, Eigen::Dynamic, 3>& v,
Eigen::Matrix<Real, Eigen::Dynamic, 3>& a) {
	for(unsigned int i = 0; i < N; i++) {
		Eigen::Vector3f x = r.row(i).reshaped(3, 1);

		Eigen::Vector3f c = G_SI * -x.normalized() * 1000.0 * m(i) * (1.0 / x.dot(x));

		a(i, 0) = c(0);
		a(i, 1) = c(1);
		a(i, 2) = c(2);
	}

	/*
	// Compute Forces (O(N^2), very slow...)
	{
		for(unsigned int i = 0; i < N; i++) {
			float m1 = m[i];

			float xi[3] = {{x[3 * i + 0]}, {x[3 * i + 1]}, {x[3 * i + 2]}};

			//#pragma omp parallel for
			for(unsigned int j = 0; j < N; j++) {
				if(i == j) continue;

				float m2 = m[j];

				float d[3] = {
					{x[3 * j + 0] - xi[0]},
					{x[3 * j + 1] - xi[1]},
					{x[3 * j + 2] - xi[2]}
				};

				float r2 = (d[0] * d[0]) + (d[1] * d[1]) + (d[2] * d[2]);
				float r  = sqrtf(r2);

				a[3 * i + 0] += d[0] * G * m1 * m2 * ((float)1 / (sqrtf(r) * r2));
				a[3 * i + 1] += d[1] * G * m1 * m2 * ((float)1 / (sqrtf(r) * r2));
				a[3 * i + 2] += d[2] * G * m1 * m2 * ((float)1 / (sqrtf(r) * r2));
			}
		}
	}
	*/
}

void Render(
unsigned int N,
Eigen::Matrix<Real, Eigen::Dynamic, 1>& m,
Eigen::Matrix<Real, Eigen::Dynamic, 3>& r,
Eigen::Matrix<Real, Eigen::Dynamic, 3>& v,
Eigen::Matrix<Real, Eigen::Dynamic, 3>& a) {
}

int main(int argc, char** argv) {
	if(argc > 1) {
		std::cerr << "Invalid argument: \"" << argv[1] << "\"" << std::endl;

		return EXIT_FAILURE;
	}

	unsigned int N = 65536; // number of particles

	//double dt = 0.06666666666666666666666666666667;
	double dt = 0.66666666666666666666666666666667;

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
	//float* m = (float*)malloc(sizeof(float) * N * 1); // Particle masses
	//float* x = (float*)malloc(sizeof(float) * N * 3); // Particle positions
	//float* v = (float*)malloc(sizeof(float) * N * 3); // Particle velocities
	//float* a = (float*)malloc(sizeof(float) * N * 3); // Particle accelerations

	Eigen::Matrix<Real, Eigen::Dynamic, 1> m = Eigen::MatrixXf::Constant(N, 1, 1.0 / (double)N);
	Eigen::Matrix<Real, Eigen::Dynamic, 3> r = Eigen::MatrixXf::Zero(N, 3);
	Eigen::Matrix<Real, Eigen::Dynamic, 3> v = Eigen::MatrixXf::Zero(N, 3);
	Eigen::Matrix<Real, Eigen::Dynamic, 3> a = Eigen::MatrixXf::Zero(N, 3);

	std::srand(1);

	for(unsigned int i = 0; i < N; i++) {
		r(i, 0) = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		r(i, 1) = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		r(i, 2) = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
	}

	std::srand(2);

	for(unsigned int i = 0; i < N; i++) {
		v(i, 0) = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		v(i, 1) = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		v(i, 2) = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
	}

	UpdateAcceleration(N, m, r, v, a); // update acceleration

	double t = 0.0;

	for(int n = 1; n < 300; n++) {
		//#pragma omp simd
		for(unsigned int i = 0; i < (ImageSizeX * ImageSizeY); i++) image.data[i] = RGBA32F(0.0, 0.0, 0.0, 1.0); // Clear image buffer

		// Render the image
		for(unsigned int i = 0; i < N; i++) {
			// vec2 uv = (fragCoord - 0.5 * iResolution.xy) / min(iResolution.x, iResolution.y);
			// vec2 fragCoord = (min(iResolution.x, iResolution.y) * uv) + (0.5 * iResolution.xy);

			Eigen::Vector3f p = r.row(i);

			double _t = 0.0125 * t;

			Eigen::Vector3f _p(
			(Real)(p(0) * cos(2.0 * acos(-1.0) * _t)) - (p(2) * sin(2.0 * acos(-1.0) * _t)),
			(Real)(p(0) * sin(2.0 * acos(-1.0) * _t)) + (p(2) * sin(2.0 * acos(-1.0) * _t)),
			p(1)
			);

			p = _p;

			float uv_x = (fminf(ImageSizeX, ImageSizeY) * (3e-2 * p(0))) + (0.5 * ImageSizeX);
			float uv_y = (fminf(ImageSizeX, ImageSizeY) * (3e-2 * p(2))) + (0.5 * ImageSizeY);

			int coord_x = (int)uv_x;
			int coord_y = (int)uv_y;

			if(0 <= coord_x && coord_x < ImageSizeX && 0 <= coord_y && coord_y < ImageSizeY) {
				image.data[coord_x + ImageSizeX * coord_y].r += 0.0625;
				image.data[coord_x + ImageSizeX * coord_y].g += 0.0625;
				image.data[coord_x + ImageSizeX * coord_y].b += 0.0625;
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
		//for(unsigned int i = 0; i < 3 * N; i++) v[i] += 0.5 * dt * a[i]; // 1/2 kick
		v = v + (0.5 * dt * a); // 1/2 kick

		//#pragma omp simd
		//for(unsigned int i = 0; i < 3 * N; i++) x[i] += dt * v[i]; // drift
		r = r + (dt * v); // drift

		UpdateAcceleration(N, m, r, v, a); // update acceleration

		//#pragma omp simd
		//for(unsigned int i = 0; i < 3 * N; i++) v[i] += 0.5 * dt * a[i]; // 1/2 kick
		v = v + (0.5 * dt * a); // 1/2 kick

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
