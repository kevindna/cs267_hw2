#include "common.h"
#include <cmath>
#include <array>
#include <vector>

#define NUM_TILES_DIM  16
#define NUM_TILES  		NUM_TILES_DIM*NUM_TILES_DIM
#define GRID_DIM  		(NUM_TILES_DIM + 2)
#define GRID_SIZE 		 GRID_DIM*GRID_DIM + 2*GRID_DIM
//#define ORG

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

// Partition all particles into tiles
void partition_particles(particle_t* p, int num_p, double size) {
	int row, col, ind;


	for(int i = 0; i < NUM_TILES; i++) {
		tiles[i].clear();
		bin_cnt[i] = 0;
	}

//#pragma GCC unroll 10
	for(int i = 0; i < num_p; i++) {
		row = floor(p[i].y/NUM_TILES_DIM);
		col = floor(p[i].x/NUM_TILES_DIM);
		ind = row*(NUM_TILES_DIM+2) + col+1;

		tiles[ind].push_back(&p[i]); // = i;
		bin_cnt[ind]++;
	}

	for(int r = 0; r < GRID_DIM; r++) {
		rows[r].clear();
		for(int c = 0; c < GRID_DIM; c++) {
			rows[r].insert(rows[r].end(), tiles[r*GRID_DIM + c].begin(), tiles[r*GRID_DIM + c].end());
		}
	}

}

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
	// that you may need. This function will be called once before the
	// algorithm begins. Do not do any particle simulation here

	// Reserve size to prevent to many dynamic allocation (memcpy calls)
	for (int x = 0; x < GRID_SIZE; x++) {
		tiles[x].reserve(num_parts/4);
	}

	// partition_particles(parts, num_parts, size);
}


static inline int get_tile_ind(int x, int y) {
	int row = floor(y/NUM_TILES_DIM);
	int col = floor(x/NUM_TILES_DIM);
	return row*GRID_DIM + GRID_DIM + col + 1;
}

/*
 * Runt he apply_force function across an entire tile
 *
 * */
static inline void apply_force_tile(particle_t* parts, particle_t p, std::vector<particle_t *> tile, int tile_num_parts) {
	for (int x = 0; x < tile_num_parts; x++) {
													//apply_force(p, parts[tile[x]]);
													apply_force(p, *tile[x]);
												}
}


// Original function
void simulate_one_step(particle_t* parts, int num_parts, double size) {
#ifdef ORG
	// Compute Forces
	for (int i = 0; i < num_parts; ++i) {
		parts[i].ax = parts[i].ay = 0;
		for (int j = 0; j < num_parts; ++j) {
			apply_force(parts[i], parts[j]);
		}
	}

	// Move Particles
	for (int i = 0; i < num_parts; ++i) {
		move(parts[i], size);
	}

#else
	const int neighbor0_off = -GRID_DIM - 1;    // Upper left
	const int neighbor1_off = -GRID_DIM;	    // Upper
	const int neighbor2_off = -GRID_DIM + 1;    // Upper right
	const int neighbor3_off = - 1; 				// Immediate left
	const int self_off = 0;
	const int neighbor5_off =  1; 			    // Immediate right
	const int neighbor6_off = GRID_DIM - 1; 	// Lower left
	const int neighbor7_off = GRID_DIM; 		// Bottom
	const int neighbor8_off = GRID_DIM + 1; 	// Bottom right

	partition_particles(parts, num_parts, size);;

	for (int i = 0; i < num_parts; ++i) {
		int tile_ind = get_tile_ind(parts[i].x, parts[i].y);
		parts[i].ax = parts[i].ay = 0; 	// Clear acceleration
		apply_force_tile(parts, parts[i], tiles[tile_ind], bin_cnt[tile_ind]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor0_off], bin_cnt[tile_ind + neighbor0_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor1_off], bin_cnt[tile_ind + neighbor1_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor2_off], bin_cnt[tile_ind + neighbor2_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor3_off], bin_cnt[tile_ind + neighbor3_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor5_off], bin_cnt[tile_ind + neighbor5_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor6_off], bin_cnt[tile_ind + neighbor6_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor7_off], bin_cnt[tile_ind + neighbor7_off]);
		apply_force_tile(parts, parts[i], tiles[tile_ind + neighbor8_off], bin_cnt[tile_ind + neighbor8_off]);
	}

	// Move Particles
	for (int i = 0; i < num_parts; ++i) {
		move(parts[i], size);
	}

#endif
}


/*
	 void add_neighbors(std::vector<int>* neighbors, int x, int y) {
	 if (x >= 0 && y >= 0) {
	 neighbors->insert(neighbors->end(), tiles[x*bin_dim + y].begin(), tiles[y*bin_dim + x].end());
	 }
	 }
 */


	// Compute Forces
	/*
		 for (int tile = 0; tile < NUM_TILES; tile++) {
		 for (int p = 0; p < bin_cnt[tile]; p++) {
		 particle_t particle_A = parts[tiles[tile][p]];
		 particle_A.x = 0;
		 particle_A.y = 0;
		 for (int q = 0; q < bin_cnt[tile]; q++) {
		 particle_t particle_B = parts[tiles[tile][q]];
		 apply_force(particle_A, particle_B);
		 }
		 }
		 }
	 */
	//  for (int i = 0; i < num_parts; ++i) {
	//		int grid_ind = get_tile_ind(parts[i].x, parts[i].y);
	//
	//	}

		//
		//				int row = floor(parts[i].y/bin_dim);
		//				int col = floor(parts[i].x/bin_dim);
		//
		//				std::vector<int> neighbors = tiles[row*bin_dim + col];
		///*				neighbors.insert(neighbors.end(), tiles[(row-1)*bin_dim + col-1].begin(), tiles[(row-1)*bin_dim + col-1].end());
		//				neighbors.insert(neighbors.end(), tiles[(row+0)*bin_dim + col-1].begin(), tiles[(row+0)*bin_dim + col-1].end());
		//				neighbors.insert(neighbors.end(), tiles[(row+1)*bin_dim + col-1].begin(), tiles[(row+1)*bin_dim + col-1].end());
		//				neighbors.insert(neighbors.end(), tiles[(row+1)*bin_dim + col+0].begin(), tiles[(row+1)*bin_dim + col+0].end());
		//				neighbors.insert(neighbors.end(), tiles[(row-1)*bin_dim + col+0].begin(), tiles[(row-1)*bin_dim + col+0].end());
		//				neighbors.insert(neighbors.end(), tiles[(row-1)*bin_dim + col+1].begin(), tiles[(row-1)*bin_dim + col+1].end());
		//				neighbors.insert(neighbors.end(), tiles[(row+0)*bin_dim + col+1].begin(), tiles[(row+0)*bin_dim + col+1].end());
		//				neighbors.insert(neighbors.end(), tiles[(row+1)*bin_dim + col+1].begin(), tiles[(row+1)*bin_dim + col+1].end());
		//*/
		//			add_neighbors(&neighbors, (row-1), col-1);
		//			add_neighbors(&neighbors, (row+0), col-1);
		//			add_neighbors(&neighbors, (row+1), col-1);
		//			add_neighbors(&neighbors, (row+1), col+0);
		//			add_neighbors(&neighbors, (row-1), col+0);
		//			add_neighbors(&neighbors, (row-1), col+1);
		//			add_neighbors(&neighbors, (row+0), col+1);
		//			add_neighbors(&neighbors, (row+1), col+1);
		//
		//int neighbor_cnt = bin_cnt[(row-1)*bin_dim + col-1] +
		//bin_cnt[(row+0)*bin_dim + col-1] +
		//bin_cnt[(row+1)*bin_dim + col-1] +
		//bin_cnt[(row+1)*bin_dim + col+0] +
		//bin_cnt[(row-1)*bin_dim + col+0] +
		//bin_cnt[(row-1)*bin_dim + col+1] +
		//bin_cnt[(row+0)*bin_dim + col+1] +
		//bin_cnt[(row+1)*bin_dim + col+1];
		//
		//        for (int j = 0; j < neighbor_cnt; ++j) {
		//            apply_force(parts[i], parts[neighbors[j]]);
		//        }
