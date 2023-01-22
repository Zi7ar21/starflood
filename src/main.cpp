#include <iostream>

#include <cstdlib>
#include <cstdint>
#include <cmath>

//#include <array>
#include <vector>

#include <omp.h>

#include "image.hpp"

//#include <cblas.h>

//#include <image.hpp>

//#define STARFLOOD_USE_DOUBLE

#define G 10.0

#ifdef STARFLOOD_USE_DOUBLE
typedef double Real;
#else
typedef float Real;
#endif

namespace Starflood {
	typedef uint32_t DataID_t;
	typedef uint32_t LeafID_t;
	typedef uint32_t NodeID_t;

	typedef uint8_t Depth_t;

	typedef float Real;
	//typedef double Real;
	//typedef long double Real;

	class Vec2 {
		public:
			Real x, y;

			Vec2(                ): x( 0), y( 0) {}
			Vec2(Real _x         ): x(_x), y(_x) {}
			Vec2(Real _x, Real _y): x(_x), y(_y) {}

			Vec2 operator*(const Real& _x) const { return Vec2(x * _x  , y * _x  ); }
			Vec2 operator*(const Vec2& _v) const { return Vec2(x * _v.x, y * _v.y); }
			Vec2 operator/(const Real& _x) const { return Vec2(x / _x  , y / _x  ); }
			Vec2 operator/(const Vec2& _v) const { return Vec2(x / _v.x, y / _v.y); }
			Vec2 operator+(const Real& _x) const { return Vec2(x + _x  , y + _x  ); }
			Vec2 operator+(const Vec2& _v) const { return Vec2(x + _v.x, y + _v.y); }
			Vec2 operator-(const Real& _x) const { return Vec2(x - _x  , y - _x  ); }
			Vec2 operator-(const Vec2& _v) const { return Vec2(x - _v.x, y - _v.y); }
	};

	class Vec3 {
		public:
			Real x, y, z;

			Vec3(                         ): x( 0), y( 0), z( 0) {}
			Vec3(Real _x                  ): x(_x), y(_x), z(_x) {}
			Vec3(Real _x, Real _y, Real _z): x(_x), y(_y), z(_z) {}

			Vec3 operator*(const Real& _x) const { return Vec3(x * _x  , y * _x  , z * _x  ); }
			Vec3 operator*(const Vec3& _v) const { return Vec3(x * _v.x, y * _v.y, z * _v.z); }
			Vec3 operator/(const Real& _x) const { return Vec3(x / _x  , y / _x  , z / _x  ); }
			Vec3 operator/(const Vec3& _v) const { return Vec3(x / _v.x, y / _v.y, z / _v.z); }
			Vec3 operator+(const Real& _x) const { return Vec3(x + _x  , y + _x  , z + _x  ); }
			Vec3 operator+(const Vec3& _v) const { return Vec3(x + _v.x, y + _v.y, z + _v.z); }
			Vec3 operator-(const Real& _x) const { return Vec3(x - _x  , y - _x  , z - _x  ); }
			Vec3 operator-(const Vec3& _v) const { return Vec3(x - _v.x, y - _v.y, z - _v.z); }
	};

	class Vec4 {
		public:
			Real x, y, z, w;

			Vec4(                                  ): x( 0), y( 0), z( 0), w( 0) {}
			Vec4(Real _x                           ): x(_x), y(_x), z(_x), w(_x) {}
			Vec4(Real _x, Real _y, Real _z, Real _w): x(_x), y(_y), z(_z), w(_w) {}

			Vec4 operator*(const Real& _x) const { return Vec4(x * _x  , y * _x  , z * _x  , w * _x  ); }
			Vec4 operator*(const Vec4& _v) const { return Vec4(x * _v.x, y * _v.y, z * _v.z, w * _v.w); }
			Vec4 operator/(const Real& _x) const { return Vec4(x / _x  , y / _x  , z / _x  , w / _x  ); }
			Vec4 operator/(const Vec4& _v) const { return Vec4(x / _v.x, y / _v.y, z / _v.z, w / _v.w); }
			Vec4 operator+(const Real& _x) const { return Vec4(x + _x  , y + _x  , z + _x  , w + _x  ); }
			Vec4 operator+(const Vec4& _v) const { return Vec4(x + _v.x, y + _v.y, z + _v.z, w + _v.w); }
			Vec4 operator-(const Real& _x) const { return Vec4(x - _x  , y - _x  , z - _x  , w - _x  ); }
			Vec4 operator-(const Vec4& _v) const { return Vec4(x - _v.x, y - _v.y, z - _v.z, w - _v.w); }
	};

	Real dot(Vec2 _a, Vec2 _b) { return (_a.x * _b.x) + (_a.y * _b.y); }
	Real dot(Vec3 _a, Vec3 _b) { return (_a.x * _b.x) + (_a.y * _b.y) + (_a.z * _b.z); }
	Real dot(Vec4 _a, Vec4 _b) { return (_a.x * _b.x) + (_a.y * _b.y) + (_a.z * _b.z) + (_a.w * _b.w); }

	Real dot2(Vec2 _v) { return (_v.x * _v.x) + (_v.y * _v.y); }
	Real dot2(Vec3 _v) { return (_v.x * _v.x) + (_v.y * _v.y) + (_v.z * _v.z); }
	Real dot2(Vec4 _v) { return (_v.x * _v.x) + (_v.y * _v.y) + (_v.z * _v.z) + (_v.w * _v.w); }

	#ifdef STARFLOOD_USE_DOUBLE
	Real length(Vec2 _v) { return (Real)sqrtf((float)dot2(_v)); }
	Real length(Vec3 _v) { return (Real)sqrtf((float)dot2(_v)); }
	Real length(Vec4 _v) { return (Real)sqrtf((float)dot2(_v)); }
	#else
	Real length(Vec2 _v) { return (Real)sqrt((double)dot2(_v)); }
	Real length(Vec3 _v) { return (Real)sqrt((double)dot2(_v)); }
	Real length(Vec4 _v) { return (Real)sqrt((double)dot2(_v)); }
	#endif

	class Body {
		private:
		public:
			Real m;
			Real h;
			Real r[3];
			Real v[3];

			Body() {
			}
	};

	class Data {
		private:
		public:
	};

	class Node {
		private:
			bool initialized;
			LeafID_t parent;
			LeafID_t children;
			//std::vector<DataID_t>;
		public:
			Real x[]; // node center of mass
	};

	class Tree {
		private:
			uint8_t depthMinimum;
			uint8_t depthMaximum;
			std::vector<Data> data;
			std::vector<Node> node;
		public:
			Tree() {
			}

			void Insert() {
				Depth_t depth;

				while(true) {
					if(depth >= depthMaximum) {
					}
				}
			}

			void Remove() {
				Depth_t depth;

				while(true) {
					if(depth >= depthMaximum) {
					}
				}
			}

			void Search() {
				Depth_t depth;

				while(true) {
				}
			}
	};

	class Mesh {
		private:
		public:
			Mesh() {
			}
	};

	class Simulation {
		private:
			uint32_t N;
			bool initialized;
		public:
			//std::vector<Vec3> x; // position
			//std::vector<Vec3> v; // velocity

			Simulation() {
				initialized = false;
			}

			void Initialize(uint32_t _N) {
				initialized = true;

				N = _N;
			}

			void Deallocate() {
				initialized = false;
			}
	};

	/*
	void UpdateSimulation(Simulation& simulation) {
		std::vector<std::array<float, 3>> x = simulation.x;
		std::vector<Vec3> v = simulation.v;

		for(int i = 0; i < simulation.N; i++) {
			Vec3 F = Vec3(0);

			x[i] = x[i] + F;
		}

		simulation.x = x;
		simulation.v = v;
	}
	*/
}

//using namespace Starflood;

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

float W(float x[3], float h) {
}

void CalculateAcceleration(size_t N, float* m, float* x, float* v, float *a) {
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
}

int main(int argc, char** argv) {
	if(argc > 1) {
		std::cerr << "Invalid argument: \"" << argv[1] << "\"" << std::endl;

		return EXIT_FAILURE;
	}

	unsigned int N = 1024; // number of particles

	double dt = 0.06666666666666666666666666666667;

	unsigned int ImageSizeX = 256;
	unsigned int ImageSizeY = 256;

	std::cout
	<< "+-----------+\n"
	<< "| Starflood |\n"
	<< "+-----------+\n"
	<< "\n"
	<< "Version 0.0.1\n"
	<< "\n"
	<< "Parameters\n"
	<< "{\n"
	<< "N: " << N << "\n"
	<< "\n"
	<< "\n"
	<< "}\n"
	<< std::endl;

	//return EXIT_SUCCESS;

	Buffer2D<RGBA32F> image(ImageSizeX, ImageSizeY);

	//simulation.Initialize();

	// row-major
	float* m = (float*)calloc(1 * N, sizeof(float)); // Particle masses
	float* x = (float*)calloc(3 * N, sizeof(float)); // Particle positions
	float* v = (float*)calloc(3 * N, sizeof(float)); // Particle velocities
	float* a = (float*)calloc(3 * N, sizeof(float)); // Particle accelerations

	std::srand(1);

	for(unsigned int i = 0; i < N; i++) {
		m[i] = 1.0 / (double)N;

		x[3 * i + 0] = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		x[3 * i + 1] = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		x[3 * i + 2] = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);

		v[3 * i + 0] = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		v[3 * i + 1] = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
		v[3 * i + 2] = 0.1 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);

		a[3 * i + 0] = 0.0;
		a[3 * i + 1] = 0.0;
		a[3 * i + 2] = 0.0;
	}

	CalculateAcceleration(N, m, x, v, a);

	/*
	//#pragma omp simd
	for(int i = 0; i < N; i++) {
		x[3 * i + 0] = (Real)0;
		x[3 * i + 1] = (Real)0;
		x[3 * i + 2] = (Real)0;

		v[3 * i + 0] = (Real)0;
		v[3 * i + 1] = (Real)0;
		v[3 * i + 2] = (Real)0;
	}
	*/

	double t = 0.0;

	for(int n = 0; n < 150; n++) {
		//#pragma omp simd
		for(unsigned int i = 0; i < (ImageSizeX * ImageSizeY); i++) image.data[i] = RGBA32F(0.0, 0.0, 0.0, 1.0); // Clear image buffer

		// Render the image
		for(unsigned int i = 0; i < N; i++) {
			// vec2 uv = (fragCoord - 0.5 * iResolution.xy) / min(iResolution.x, iResolution.y);
			// vec2 fragCoord = (min(iResolution.x, iResolution.y) * uv) + (0.5 * iResolution.xy);

			float  p[3] = {x[3 * i + 0], x[3 * i + 1], x[3 * i + 2]};

			double _t = 0.05 * t;

			float _p[3] = {
			(Real)(p[0] * cos(2.0 * acos(-1.0) * _t)) - (p[2] * sin(2.0 * acos(-1.0) * _t)),
			(Real)(p[0] * sin(2.0 * acos(-1.0) * _t)) + (p[2] * sin(2.0 * acos(-1.0) * _t)),
			p[1]
			};

			p[0] = _p[0];
			p[1] = _p[2];
			p[2] = _p[1];

			float uv_x = (fminf(ImageSizeX, ImageSizeY) * (2e-2 * p[0])) + (0.5 * ImageSizeX);
			float uv_y = (fminf(ImageSizeX, ImageSizeY) * (2e-2 * p[1])) + (0.5 * ImageSizeY);

			int coord_x = (int)uv_x;
			int coord_y = (int)uv_y;

			if(0 <= coord_x && coord_x < ImageSizeX && 0 <= coord_y && coord_y < ImageSizeY) {
				image.data[coord_x + ImageSizeX * coord_y].r += 0.5;
				image.data[coord_x + ImageSizeX * coord_y].g += 0.5;
				image.data[coord_x + ImageSizeX * coord_y].b += 0.5;
			}
		}

		// Write the image buffer
		{
			char filename[80];

			sprintf(filename, "out/IMG_%04d.png", n);
			//sprintf(filename, "out/IMG_%04d.hdr", n);

			WriteImageBuffer(image, filename);
			//WriteImageBufferHDR(image, filename);
		}

		//UpdateSimulation(simulation);

		/*
		#pragma omp simd
		for(unsigned int i = 0; i < N; i++) {
			// 1/2 kick
			v[3 * i + 0] += 0.5 * dt * a[3 * i + 0];
			v[3 * i + 1] += 0.5 * dt * a[3 * i + 1];
			v[3 * i + 2] += 0.5 * dt * a[3 * i + 2];
		}
		*/

		//#pragma omp simd
		for(unsigned int i = 0; i < 3 * N; i++) v[i] += 0.5 * dt * a[i]; // 1/2 kick

		/*
		#pragma omp simd
		for(unsigned int i = 0; i < N; i++) {
			// drift
			x[3 * i + 0] += dt * v[3 * i + 0];
			x[3 * i + 1] += dt * v[3 * i + 1];
			x[3 * i + 2] += dt * v[3 * i + 2];
		}
		*/

		//#pragma omp simd
		for(unsigned int i = 0; i < 3 * N; i++) x[i] += dt * v[i]; // drift

		CalculateAcceleration(N, m, x, v, a);

		/*
		#pragma omp simd
		for(unsigned int i = 0; i < N; i++) {
			// 1/2 kick
			v[3 * i + 0] += 0.5 * dt * a[3 * i + 0];
			v[3 * i + 1] += 0.5 * dt * a[3 * i + 1];
			v[3 * i + 2] += 0.5 * dt * a[3 * i + 2];
		}
		*/

		//#pragma omp simd
		for(unsigned int i = 0; i < 3 * N; i++) v[i] += 0.5 * dt * a[i]; // 1/2 kick

		t += dt;
	}

	free(m);
	free(x);
	free(v);
	free(a);

	//simulation.Deallocate();

	std::cout << "Finished!" << std::endl;

	return EXIT_SUCCESS;
}
