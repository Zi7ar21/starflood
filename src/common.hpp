#pragma once

#include "config.h"

#ifdef STARFLOOD_USE_DOUBLE
typedef double Real;
#else
typedef float Real;
#endif

namespace Starflood {
	//typedef uint32_t DataID_t;
	//typedef uint32_t LeafID_t;
	//typedef uint32_t NodeID_t;

	//typedef uint8_t Depth_t;

	/*
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
			uint32_t id;
			Real m;
			Real h;
			Real r[3];
			Real v[3];
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
			//std::vector<Data> data;
			//std::vector<Node> node;
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
			bool initialized;
		public:
			uint32_t N;
			//std::vector<Vec3> x; // position
			//std::vector<Vec3> v; // velocity
			//std::vector<std::array<Real, 3>> x;

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
	*/

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

	class Tree {
		private:
			//uint32_t ;

		public:
			Tree() {
			}
	};
}
