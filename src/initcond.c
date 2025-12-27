#include "initcond.h"

#include <math.h>

#include "common.h"
#include "config.h"
#include "rng.h"
#include "types.h"

int initcond_generate(real* restrict mas, real* restrict rad, real* restrict pos, real* restrict vel, unsigned int N) {
	const double body_mass = 1.0 / (double)N;

	// Initialize Mass
	for(unsigned int i = 0u; i < N; i++) {
		mas[i] = (real)body_mass;
	}

	#ifdef ENABLE_SPH
	// Initialize Smoothing Radius
	for(unsigned int i = 0u; i < N; i++) {
		rad[i] = (real)0.100;
	}
	#endif

	// Initialize Position
	for(unsigned int i = 0u; i < N; i++) {
		double p[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		uint32_t s[4] = {
			(uint32_t)0xD88D5EFDu + (uint32_t)i,
			(uint32_t)0x203F83D3u,
			(uint32_t)0x710BE493u,
			(uint32_t)0xA294E7F1u
		};

		pcg4d(s);
		//pcg4d(s); // second round for better statistical quality

		{
			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			double the = acos(2.0*r[0]-1.0), phi = TAU * r[1];

			double cos_the = cos(the), sin_the = sin(the);
			double cos_phi = cos(phi), sin_phi = sin(phi);

			// random point on unit sphere
			p[0u] = cos_phi * sin_the;
			p[1u] = sin_phi * sin_the;
			p[2u] =           cos_the;

			// random point in unit ball
			p[0] *= cbrt(r[2]);
			p[1] *= cbrt(r[2]);
			p[2] *= cbrt(r[2]);
		}

		{
			pcg4d(s);

			double r[4] = {
				INV_PCG32_MAX * (double)s[0],
				INV_PCG32_MAX * (double)s[1],
				INV_PCG32_MAX * (double)s[2],
				INV_PCG32_MAX * (double)s[3]
			};

			// Standard Normal Distribution
			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			double n[4] = {
				sqrt( -2.0 * log(r[0]) ) * cos(TAU * r[1]),
				sqrt( -2.0 * log(r[0]) ) * sin(TAU * r[1]),
				sqrt( -2.0 * log(r[2]) ) * cos(TAU * r[3]),
				sqrt( -2.0 * log(r[2]) ) * sin(TAU * r[3])
			};

			//p[0] = 1.000 * n[0];
			p[1] = 1.000 * n[1];
			//p[2] = 1.000 * n[2];
		}

		/*
		p[0] *= 1.000 * sqrt(r[2]);
		p[1] *= 0.025 * sqrt(r[2]);
		p[2] *= 1.000 * sqrt(r[2]);
		*/

		p[0] *= 1.000;
		p[1] *= 0.025;
		p[2] *= 1.000;

		pos[3u*i+0u] = (real)p[0u];
		pos[3u*i+1u] = (real)p[1u];
		pos[3u*i+2u] = (real)p[2u];
	}

	// Initialize Velocity
	for(unsigned int i = 0u; i < N; i++) {
		double p[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		double v[3] = {
			vel[3u*i+0u],
			vel[3u*i+1u],
			vel[3u*i+2u]
		};

		uint32_t s[4] = {
			(uint32_t)0x4AFDC7CFu + (uint32_t)i,
			(uint32_t)0xAE8D7CF2u,
			(uint32_t)0xE4827F00u,
			(uint32_t)0xE44CA389u
		};

		pcg4d(s);
		//pcg4d(s); // second round for better statistical quality

		double r[4] = {
			INV_PCG32_MAX * (double)s[0],
			INV_PCG32_MAX * (double)s[1],
			INV_PCG32_MAX * (double)s[2],
			INV_PCG32_MAX * (double)s[3]
		};

		// Standard Normal Distribution
		// Box-Muller Transform
		// https://en.wikipedia.org/wiki/Box–Muller_transform
		double n[4] = {
			sqrt( -2.0 * log(r[0]) ) * cos(TAU * r[1]),
			sqrt( -2.0 * log(r[0]) ) * sin(TAU * r[1]),
			sqrt( -2.0 * log(r[2]) ) * cos(TAU * r[3]),
			sqrt( -2.0 * log(r[2]) ) * sin(TAU * r[3])
		};

		/*
		double r2 = (p[0]*p[0])+(p[2]*p[2]);

		double inv_r2 = 1.0 / (      r2+0.001  );
		double inv_r1 = 1.0 / ( sqrt(r2+0.001) );

		
		v[0] = 0.250 * sqrt(G * body_mass) *  p[2] + 1.000e-6 * n[0];
		v[1] = 0.000 * sqrt(G * body_mass) *  p[1] + 1.000e-6 * n[1];
		v[2] = 0.250 * sqrt(G * body_mass) * -p[0] + 1.000e-6 * n[2];
		*/


		v[0] = 5.000e-3 * sqrt(G * 1.0) *  p[2] + 1.000e-6 * n[0];
		v[1] = 0.000e-3 * sqrt(G * 1.0) *  p[1] + 8.000e-6 * n[1];
		v[2] = 5.000e-3 * sqrt(G * 1.0) * -p[0] + 1.000e-6 * n[2];

		/*
		v[0] = 3.000e-3 * p[0] + 1.000e-6 * n[0];
		v[1] = 3.000e-3 * p[1] + 1.000e-6 * n[1];
		v[2] = 3.000e-3 * p[2] + 1.000e-6 * n[2];
		*/

		/*
		v[0] = 1.0e-6 * n[0];
		v[1] = 1.0e-6 * n[1];
		v[2] = 1.0e-6 * n[2];
		*/

		vel[3u*i+0u] = (real)v[0u];
		vel[3u*i+1u] = (real)v[1u];
		vel[3u*i+2u] = (real)v[2u];
	}

	// Transform Position/Velocity
	for(unsigned int i = 0u; i < N; i++) {
		double p[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		double v[3] = {
			vel[3u*i+0u],
			vel[3u*i+1u],
			vel[3u*i+2u]
		};

		uint32_t s[4] = {
			(uint32_t)0x5D20CDF7u + (uint32_t)i,
			(uint32_t)0xE5C4E174u,
			(uint32_t)0xCB77B51Eu,
			(uint32_t)0xA82C87A3u
		};

		pcg4d(s);
		//pcg4d(s); // second round for better statistical quality

		double r[4] = {
			INV_PCG32_MAX * (double)s[0],
			INV_PCG32_MAX * (double)s[1],
			INV_PCG32_MAX * (double)s[2],
			INV_PCG32_MAX * (double)s[3]
		};

		// Thin disk
		//if(0.000 <= r[0] && r[0] < 0.850) {
		//	p[1] *= 0.025;
		//}

		//p[0] += r[1] < 0.500 ?  4.000 : -4.000;
		/*
		p[1] += r[1] < 0.500 ?  1.000 : -1.000;
		v[0] *= r[1] < 0.500 ?  1.000 : -1.000;
		v[2] *= r[1] < 0.500 ?  1.000 : -1.000;
		v[0] += r[1] < 0.500 ? -0.005 :  0.005;
		*/
		//v[1] += r[1] < 0.500 ? -0.005 :  0.005;

		pos[3u*i+0u] = (real)p[0u];
		pos[3u*i+1u] = (real)p[1u];
		pos[3u*i+2u] = (real)p[2u];

		vel[3u*i+0u] = (real)v[0u];
		vel[3u*i+1u] = (real)v[1u];
		vel[3u*i+2u] = (real)v[2u];
	}

	return STARFLOOD_SUCCESS;
}
