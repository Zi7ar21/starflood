#include "initcond.h"

#include <math.h>

#include "common.h"
#include "config.h"
#include "rng.h"
#include "types.h"

int initcond_generate(real* volatile mas, real* volatile pos, real* volatile vel, unsigned int N) {
	const double body_mass = 1.0 / (double)N;

	// Initialize Mass
	for(unsigned int i = 0u; i < N; i++) {
		mas[i] = (real)body_mass;
	}

	// Initialize Position
	for(unsigned int i = 0u; i < N; i++) {
		real p[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		uint32_t s[4] = {
			(uint32_t)i,
			(uint32_t)0x203F83D3u,
			(uint32_t)0x710BE493u,
			(uint32_t)0xA294E7F1u
		};

		pcg4d(s);
		pcg4d(s); // second round for better statistical quality

		double r[4] = {
			INV_PCG32_MAX * (double)s[0],
			INV_PCG32_MAX * (double)s[1],
			INV_PCG32_MAX * (double)s[2],
			INV_PCG32_MAX * (double)s[3]
		};

		// Box-Muller Transform
		// https://en.wikipedia.org/wiki/Box–Muller_transform
		double n[3] = {
			sqrt( -2.0 * log(r[0u]) ) * cos(TAU * r[2u]),
			sqrt( -2.0 * log(r[0u]) ) * sin(TAU * r[2u]),
			sqrt( -2.0 * log(r[1u]) ) * cos(TAU * r[3u])
		};

		p[0u] = 1.000 * n[0u];
		p[1u] = 1.000 * n[1u];
		p[2u] = 1.000 * n[2u];

		pos[3u*i+0u] = p[0u];
		pos[3u*i+1u] = p[1u];
		pos[3u*i+2u] = p[2u];
	}

	// Initialize Velocity
	for(unsigned int i = 0u; i < N; i++) {
		real p[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		real v[3] = {
			vel[3u*i+0u],
			vel[3u*i+1u],
			vel[3u*i+2u]
		};

		uint32_t s[4] = {
			(uint32_t)i,
			(uint32_t)0xAE8D7CF2u,
			(uint32_t)0xE4827F00u,
			(uint32_t)0xE44CA389u
		};

		pcg4d(s);
		pcg4d(s); // second round for better statistical quality

		double r[4] = {
			INV_PCG32_MAX * (double)s[0],
			INV_PCG32_MAX * (double)s[1],
			INV_PCG32_MAX * (double)s[2],
			INV_PCG32_MAX * (double)s[3]
		};

		// Box-Muller Transform
		// https://en.wikipedia.org/wiki/Box–Muller_transform
		double n[3] = {
			sqrt( -2.0 * log(r[0u]) ) * cos(TAU * r[2u]),
			sqrt( -2.0 * log(r[0u]) ) * sin(TAU * r[2u]),
			sqrt( -2.0 * log(r[1u]) ) * cos(TAU * r[3u])
		};

		/*
		double r2 = (p[0u]*p[0u])+(p[2u]*p[2u]);

		double inv_r2 = 1.0 / (      r2+0.1  );
		double inv_r1 = 1.0 / ( sqrt(r2+0.1) );
		*/

		v[0u] = 0.250 * sqrt(G * body_mass) *  p[2u] + 0.000001 * n[0u];
		v[1u] = 0.000 * sqrt(G * body_mass) *  p[1u] + 0.000001 * n[1u];
		v[2u] = 0.250 * sqrt(G * body_mass) * -p[0u] + 0.000001 * n[2u];

		vel[3u*i+0u] = v[0u];
		vel[3u*i+1u] = v[1u];
		vel[3u*i+2u] = v[2u];
	}

	// Transform Position/Velocity
	for(unsigned int i = 0u; i < N; i++) {
		real p[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		real v[3] = {
			vel[3u*i+0u],
			vel[3u*i+1u],
			vel[3u*i+2u]
		};

		uint32_t s[4] = {
			(uint32_t)i,
			(uint32_t)0xE5C4E174u,
			(uint32_t)0xCB77B51Eu,
			(uint32_t)0xA82C87A3u
		};

		pcg4d(s);
		pcg4d(s); // second round for better statistical quality

		double r[4] = {
			INV_PCG32_MAX * (double)s[0],
			INV_PCG32_MAX * (double)s[1],
			INV_PCG32_MAX * (double)s[2],
			INV_PCG32_MAX * (double)s[3]
		};

		// Thin disk
		if(0.000 <= r[0u] && r[0u] < 0.850) {
			p[1u] *= 0.025;
		}

		/*
		p[0u] += r[1u] < 0.500 ?  4.000 : -4.000;
		p[1u] += r[1u] < 0.500 ?  1.000 : -1.000;
		v[0u] *= r[1u] < 0.500 ?  1.000 : -1.000;
		v[2u] *= r[1u] < 0.500 ?  1.000 : -1.000;
		v[0u] += r[1u] < 0.500 ? -0.005 :  0.005;
		v[1u] += r[1u] < 0.500 ? -0.005 :  0.005;
		*/

		pos[3u*i+0u] = p[0u];
		pos[3u*i+1u] = p[1u];
		pos[3u*i+2u] = p[2u];

		vel[3u*i+0u] = v[0u];
		vel[3u*i+1u] = v[1u];
		vel[3u*i+2u] = v[2u];
	}

	return STARFLOOD_SUCCESS;
}
