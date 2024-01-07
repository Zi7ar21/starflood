#pragma once

#include <common.h>

#include <stack>
#include <vector>

class Node {
	public:
		/*
		Quadrants:
		+---+---+
		| 2 | 3 |
		+---+---+
		| 0 | 1 |
		+---+---+
		*/

		int parent; // index of this node's parent, -1: this is the root node
		int child0; // index of this node's child in quadrant 0
		int child1; // index of this node's child in quadrant 1
		int child2; // index of this node's child in quadrant 2
		int child3; // index of this node's child in quadrant 3

		int id; // id of the particle contained
		int np; // number of particles contained in this node and children

		real x_min; // this node's axis-aligned bounding-box x-minimum
		real y_min; // this node's axis-aligned bounding-box y-minimum
		real x_max; // this node's axis-aligned bounding-box x-maximum
		real y_max; // this node's axis-aligned bounding-box y-maximum

		real m; // mass contained in this node and children
		real x; // this node's center of mass in the x-axis
		real y; // this node's center of mass in the y-axis

		Node() {}

		Node(int _parent, real _x_min, real _y_min, real _x_max, real _y_max):
		parent(_parent), child0(-1), child1(-1), child2(-1), child3(-1), id(-1), np(0),
		x_min(_x_min), y_min(_y_min), x_max(_x_max), y_max(_y_max), m(0), x(0), y(0) {}
};

void BarnesHut(std::vector<Node> &tree, float* image, int w, int h, int* ids, real* mas, real* pos, real* acc, int N, int step_num);
