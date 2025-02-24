#include "common.h"
#include <cmath>
#include <array>
#include <vector>
#include <iostream>

#define NUM_TILES_DIM  16
#define NUM_TILES  		NUM_TILES_DIM*NUM_TILES_DIM
#define GRID_DIM  		(NUM_TILES_DIM + 2)
#define GRID_SIZE 		GRID_DIM*GRID_DIM

// Define all offesets for accessing 8 NN
#define NEIGHBOR0_OFF -GRID_DIM - 1  // Upper left
#define NEIGHBOR1_OFF -GRID_DIM	     // Upper
#define NEIGHBOR2_OFF -GRID_DIM + 1  // Upper right
#define NEIGHBOR3_OFF -1 						 // Immediate left
#define NEIGHBOR5_OFF  1 			    	 // Immediate right
#define NEIGHBOR6_OFF GRID_DIM - 1   // Lower left
#define NEIGHBOR7_OFF GRID_DIM 	 		 // Bottom
#define NEIGHBOR8_OFF GRID_DIM + 1   // Bottom right


double tile_size = 0;
int bin_cnt[GRID_SIZE] = {0};
std::vector<particle_t *> tiles[GRID_SIZE];
std::vector<particle_t *> rows[GRID_DIM];

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
	// Calculate Distance
	double dx = neighbor.x - particle.x;
	double dy = neighbor.y - particle.y;
	double r2 = dx * dx + dy * dy;

	// Check if the two particles should interact
	if (r2 > cutoff * cutoff)
		return;

	r2 = fmax(r2, min_r * min_r);
	double r = sqrt(r2);

	// Very simple short-range repulsive force
	double coef = (1 - cutoff / r) / r2 / mass;
	particle.ax += coef * dx;
	particle.ay += coef * dy;
}


// Integrate the ODE
void move(particle_t& p, double size) {
	// Slightly simplified Velocity Verlet integration
	// Conserves energy better than explicit Euler method
	p.vx += p.ax * dt;
	p.vy += p.ay * dt;
	p.x += p.vx * dt;
	p.y += p.vy * dt;

	// Bounce from walls
	while (p.x < 0 || p.x > size) {
		p.x = p.x < 0 ? -p.x : 2 * size - p.x;
		p.vx = -p.vx;
	}

	while (p.y < 0 || p.y > size) {
		p.y = p.y < 0 ? -p.y : 2 * size - p.y;
		p.vy = -p.vy;
	}
}


void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
	// that you may need. This function will be called once before the
	// algorithm begins. Do not do any particle simulation here
#if ORG == 1
		printf("Orignal Serial Implementation\n");
#else
		printf("Tiled Serial Implementation\n");
#endif

	tile_size = size/NUM_TILES_DIM;

	// Reserve size to prevent to many dynamic allocation (memcpy calls)
	for (int x = 0; x < GRID_SIZE; x++) {
		tiles[x].reserve(num_parts/25);
	}

}


static inline int get_tile_ind(double x, double y) {
	int row = floor(y/tile_size);
	int col = floor(x/tile_size);
	return  (row+1)*(GRID_DIM) + col+1;
}


// Partition all particles into tiles
void partition_particles(particle_t* p, int num_p, double size) {
	int row, col, ind;

	for(int i = 0; i < GRID_SIZE; i++) {
		tiles[i].clear();
	}

	for(int i = 0; i < num_p; i++) {
		ind = get_tile_ind(p[i].x, p[i].y);
		tiles[ind].push_back(p + i);
	}

//printf("Tile particle count: %d (%d)\n", ind, bin_cnt[ind]);

}

/*
 * Runt he apply_force function across an entire tile
 *
 * */
static inline void apply_force_tile(particle_t& p, std::vector<particle_t *>& tile) {
	for (int x = 0; x < tile.size(); x++) {
		apply_force(p, *tile[x]);
	}
}


// Original function
void simulate_one_step(particle_t* parts, int num_parts, double size) {
#if ORG==1
	// Compute Forces
	for (int i = 0; i < num_parts; ++i) {
		parts[i].ax = parts[i].ay = 0;
		for (int j = 0; j < num_parts; ++j) {
			apply_force(parts[i], parts[j]);
		}
	}

#else

	int tile_ind;
	partition_particles(parts, num_parts, size);

	for (int i = 0; i < num_parts; ++i) {
		parts[i].ax = parts[i].ay = 0; 	// Clear acceleration
		tile_ind = get_tile_ind(parts[i].x, parts[i].y);

		apply_force_tile(parts[i], tiles[tile_ind]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR0_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR1_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR2_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR3_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR5_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR6_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR7_OFF]);
		apply_force_tile(parts[i], tiles[tile_ind + NEIGHBOR8_OFF]);
	}

#endif

	// Move Particles
	for (int i = 0; i < num_parts; ++i) {
		move(parts[i], size);
	}
}
