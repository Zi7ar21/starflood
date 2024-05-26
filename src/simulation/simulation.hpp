#pragma once

#include <omp.h>

#include <cstdio>
#include <gravity/barnes-hut.hpp>
#include <simulation/particles.hpp>

class sim {
	public:
		double start_time;

		// profiling timestamps
		double t0, t1;

		size_t w; // image width
		size_t h; // image height
		size_t N; // number of bodies in the simulation

		size_t image_size; // Image Buffer Size, RGBA32F
		float* image; // Image Buffer, RGBA32F

		size_t mem_size; // simulation total buffer size
		void* mem;

		particles sim_particles;

		FILE* diagfile;

		FILE* statistics_file;

		std::vector<Node> tree;

		sim(const double &_start_time): start_time(_start_time) {
		}

		void setup();

		void init();

		void run();

		void clean();
};
