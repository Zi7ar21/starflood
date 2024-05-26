#pragma once

#include <utilities/common.h>

class particles {
	public:
		size_t mas_size; // mass buffer size
		size_t pos_size; // position buffer size
		size_t vel_size; // velocity buffer size
		size_t acc_size; // acceleration buffer size
		size_t pen_size; // potential energy buffer size
		size_t ids_size; // id buffer, contains ids of node paticles are in

		size_t particle_size;

		particles(): particle_size((11*sizeof(real))+(1*sizeof(int))) {}

		real* mas; // Mass Buffer
		real* pos; // Position Buffer
		real* vel; // Velocity Buffer
		real* acc; // Acceleration Buffer
		real* pen; // Potential Energy Buffer
		int*  ids; // id buffer, contains ids of node paticles are in
};
