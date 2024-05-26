#include <simulation/simulation.hpp>

#include <cstdio>
#include <utilities/rng.hpp>
#include <cmath>

void sim::init() {
	printf("[%.09f] Setting up initial conditions.\n", omp_get_wtime()-start_time);

	auto mas = sim_particles.mas;
	auto pos = sim_particles.pos;
	auto vel = sim_particles.vel;
	auto acc = sim_particles.acc;
	auto pen = sim_particles.pen;
	auto ids = sim_particles.ids;

	// Initialize Simulation
	{
		// random number generator state
		rng_state_t rng_state;

		// intialize random number generator
		rng_state.z0 = (uint32_t)1u;
		rng_state.z1 = (uint32_t)2u;
		rng_state.z2 = (uint32_t)3u;
		rng_state.z3 = (uint32_t)4u;

		for(size_t i = 0; i < (N * 1); i++) mas[i] = (real)((100.0*(4.0/3.0)*PI)/(double)N);
		for(size_t i = 0; i < (N * 3); i++) pos[i] = (real)0;
		for(size_t i = 0; i < (N * 3); i++) vel[i] = (real)0;
		for(size_t i = 0; i < (N * 3); i++) acc[i] = (real)0;
		for(size_t i = 0; i < (N * 1); i++) pen[i] = (real)0;

		// randomly selected positions and velocities
		for(size_t i = 0; i < N; i++) {
			// intialize random number generator
			rng_state.z0 = (uint32_t)i;
			rng_state.z1 = (uint32_t)1u;
			rng_state.z2 = (uint32_t)2u;
			rng_state.z3 = (uint32_t)3u;

			update_rng(&rng_state);

			float z0 = inv_rng32_max*rng_state.z0;
			float z1 = inv_rng32_max*rng_state.z1;
			float z2 = inv_rng32_max*rng_state.z2;
			float z3 = inv_rng32_max*rng_state.z3;

			//pos[3*i+0] = (real)((double)z0-0.5);
			//pos[3*i+1] = (real)((double)z1-0.5);
			//pos[3*i+2] = (real)((double)z2-0.5);

			if(i < N/2) {
			while((pos[3*i+0]*pos[3*i+0])+(pos[3*i+2]*pos[3*i+2]) < 0.1 * 0.1) {
				update_rng(&rng_state);

				z0 = inv_rng32_max*rng_state.z0;
				z1 = inv_rng32_max*rng_state.z1;
				z2 = inv_rng32_max*rng_state.z2;
				z3 = inv_rng32_max*rng_state.z3;

				double ang_th = (double)TAU * (double)z0;
				double ang_ph = acos(2.0 * (double)z1 - 1.0);
				double ang_rh = cbrt((double)z2);
				double sin_th = sin(ang_th);
				double cos_th = cos(ang_th);
				double sin_ph = sin(ang_ph);
				double cos_ph = cos(ang_ph);
				pos[3*i+0] = (real)(ang_rh*sin_ph*cos_th);
				pos[3*i+1] = (real)0.1*(ang_rh*sin_ph*sin_th);
				pos[3*i+2] = (real)(ang_rh*cos_ph);
			}
			} else {
				// normal distribution
				pos[3*i+0] = 0.250*sqrt(-2.0*log(z0))*cos(TAU*z2);
				pos[3*i+1] = 0.250*sqrt(-2.0*log(z0))*sin(TAU*z2);
				pos[3*i+2] = 0.250*sqrt(-2.0*log(z1))*cos(TAU*z3);
			}

			//real radasdf = powf((pos[3*i+0]*pos[3*i+0])+(pos[3*i+2]*pos[3*i+2]),0.25f);
			//pos[3*i+0] /= radasdf;
			//pos[3*i+2] /= radasdf;
			//pos[3*i+0] *= 0.200f;
			//pos[3*i+1] *= 0.005f;
			//pos[3*i+2] *= 0.200f;
			//pos[2*i+0] = 0.2*sqrt(z1)*cos(TAU*z0);
			//pos[2*i+1] = 0.2*sqrt(z1)*sin(TAU*z0);
		}

		/*
		for(size_t i = 0; i < N/2; i++) {
			uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			double ang_th = (double)TAU * (double)z0;
			double ang_ph = acos(2.0 * (double)z1 - 1.0);
			double ang_rh = cbrt((double)z2);
			double sin_th = sin(ang_th);
			double cos_th = cos(ang_th);
			double sin_ph = sin(ang_ph);
			double cos_ph = cos(ang_ph);
			pos[3*i+0] = (real)(ang_rh*sin_ph*cos_th)-(real)4.0;
			pos[3*i+1] = (real)(ang_rh*sin_ph*sin_th);
			pos[3*i+2] = (real)(ang_rh*cos_ph);
		}

		for(size_t i = N/2; i < N; i++) {
			uint32_t ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			double ang_th = (double)TAU * (double)z0;
			double ang_ph = acos(2.0 * (double)z1 - 1.0);
			double ang_rh = cbrt((double)z2);
			double sin_th = sin(ang_th);
			double cos_th = cos(ang_th);
			double sin_ph = sin(ang_ph);
			double cos_ph = cos(ang_ph);
			pos[3*i+0] = (real)(ang_rh*sin_ph*cos_th)+(real)4.0;
			pos[3*i+1] = (real)(ang_rh*sin_ph*sin_th);
			pos[3*i+2] = (real)(ang_rh*cos_ph);
		}*/

		for(size_t i = 0; i < N; i++) {
			// intialize random number generator
			rng_state.z0 = (uint32_t)i;
			rng_state.z1 = (uint32_t)4u;
			rng_state.z2 = (uint32_t)5u;
			rng_state.z3 = (uint32_t)6u;

			update_rng(&rng_state);

			float z0 = inv_rng32_max*rng_state.z0;
			float z1 = inv_rng32_max*rng_state.z1;
			float z2 = inv_rng32_max*rng_state.z2;
			float z3 = inv_rng32_max*rng_state.z3;

			// normal distribution
			vel[3*i+0] = 0.250*sqrt(-2.0*log(z0))*cos(TAU*z2);
			vel[3*i+1] = 0.250*sqrt(-2.0*log(z0))*sin(TAU*z2);
			vel[3*i+2] = 0.250*sqrt(-2.0*log(z1))*cos(TAU*z3);

			/*
			//if(i < N/2) {
				vel[3*i+0] = (real)(-PI*pos[3*i+2]);
				vel[3*i+1] = (real)(0.1*pos[3*i+1]);
				vel[3*i+2] = (real)( PI*pos[3*i+0]);

				float tempr = sqrt((vel[3*i+0]*vel[3*i+0])+(vel[3*i+2]*vel[3*i+2]));

				vel[3*i+0] /= 4.0*tempr;
				vel[3*i+2] /= 4.0*tempr;
			*/
			//} else {
			//	vel[3*i+0] = (real)(0.125*(double)pos[3*i+0]);
			//	vel[3*i+1] = (real)(0.125*(double)pos[3*i+1]);
			//	vel[3*i+2] = (real)(0.125*(double)pos[3*i+2]);
			//}
			//vel[3*i+0] = (real)0;
			//vel[3*i+1] = (real)0;
			//vel[3*i+2] = (real)0;
		}

		/*
		for(size_t i = 0; i < N/2; i++) {
			//vel[3*i+0] = (real)(-0.2*(double)pos[3*i+2]);
			//vel[3*i+1] = (real)( 0.2*(double)pos[3*i+1]);
			//vel[3*i+2] = (real)( 0.2*(double)pos[3*i+0]);
			vel[3*i+0] = (real)0.0;
			vel[3*i+1] = (real)0.0;
			//vel[3*i+2] = (real)0;
			vel[3*i+2] = (real)0;
		}
		for(size_t i = N/2; i < N; i++) {
			//vel[3*i+0] = (real)(-0.2*(double)pos[3*i+2]);
			//vel[3*i+1] = (real)( 0.2*(double)pos[3*i+1]);
			//vel[3*i+2] = (real)( 0.2*(double)pos[3*i+0]);
			vel[3*i+0] = (real)-0.0;
			vel[3*i+1] = (real)0.0;
			//vel[3*i+2] = (real)0;
			vel[3*i+2] = (real)0;
		}
		*/

		// calculate initial gravitational accelerations
		BarnesHut(tree, image, w, h, ids, mas, pos, acc, pen, N, 0);
	}

	#ifdef KEEP_CENTERED
	// center
	{
		long double avg_pos[3] = {0.0l, 0.0l, 0.0l};

		for(size_t i = 0; i < N; i++) {
			avg_pos[0] += (long double)pos[3*i+0];
			avg_pos[1] += (long double)pos[3*i+1];
			avg_pos[2] += (long double)pos[3*i+2];
		}

		avg_pos[0] /= (long double)N;
		avg_pos[1] /= (long double)N;
		avg_pos[2] /= (long double)N;

		for(size_t i = 0; i < N; i++) {
			pos[3*i+0] -= avg_pos[0];
			pos[3*i+1] -= avg_pos[1];
			pos[3*i+2] -= avg_pos[2];
		}
	}
	#endif
}
