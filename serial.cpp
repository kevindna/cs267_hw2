#include "common.h"
#include <cmath>
#include <array>
#include <vector>

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

const int NUM_TILES = 256;
std::array<std::vector<particle_t>, NUM_TILES> particle_bins;
int bin_dim;
int bin_cnt[NUM_TILES] = {0};

void partition_particles(particle_t* p, int num_p, double size) {
	int row, col;
	bin_dim = (int) ceil(size/(cutoff*10));
	int num_tiles_dim = sqrt(NUM_TILES);

	for(int i = 0; i < NUM_TILES; i++) {
		particle_bins[i].clear();
		bin_cnt[i] = 0;
	}

	for(int i = 0; i < num_p; i++) {
		row = floor(p[i].y/bin_dim);
		col = floor(p[i].x/bin_dim);
		particle_bins[row*num_tiles_dim + col].push_back(p[i]);
		bin_cnt[row*num_tiles_dim + col]++;

	}

}

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
	partition_particles(parts, num_parts, size);
}

void add_neighbors(std::vector<particle_t>* neighbors, int x, int y) {
	if (x >= 0 && y >= 0) {
		neighbors->insert(neighbors->end(), particle_bins[x*bin_dim + y].begin(), particle_bins[y*bin_dim + x].end());
	}
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
		partition_particles(parts, num_parts, size);;

		// Compute Forces
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;

				int row = floor(parts[i].y/bin_dim);
				int col = floor(parts[i].x/bin_dim);

				std::vector<particle_t> neighbors = particle_bins[row*bin_dim + col];
/*				neighbors.insert(neighbors.end(), particle_bins[(row-1)*bin_dim + col-1].begin(), particle_bins[(row-1)*bin_dim + col-1].end());
				neighbors.insert(neighbors.end(), particle_bins[(row+0)*bin_dim + col-1].begin(), particle_bins[(row+0)*bin_dim + col-1].end());
				neighbors.insert(neighbors.end(), particle_bins[(row+1)*bin_dim + col-1].begin(), particle_bins[(row+1)*bin_dim + col-1].end());
				neighbors.insert(neighbors.end(), particle_bins[(row+1)*bin_dim + col+0].begin(), particle_bins[(row+1)*bin_dim + col+0].end());
				neighbors.insert(neighbors.end(), particle_bins[(row-1)*bin_dim + col+0].begin(), particle_bins[(row-1)*bin_dim + col+0].end());
				neighbors.insert(neighbors.end(), particle_bins[(row-1)*bin_dim + col+1].begin(), particle_bins[(row-1)*bin_dim + col+1].end());
				neighbors.insert(neighbors.end(), particle_bins[(row+0)*bin_dim + col+1].begin(), particle_bins[(row+0)*bin_dim + col+1].end());
				neighbors.insert(neighbors.end(), particle_bins[(row+1)*bin_dim + col+1].begin(), particle_bins[(row+1)*bin_dim + col+1].end());
*/
			add_neighbors(&neighbors, (row-1), col-1);
			add_neighbors(&neighbors, (row+0), col-1);
			add_neighbors(&neighbors, (row+1), col-1);
			add_neighbors(&neighbors, (row+1), col+0);
			add_neighbors(&neighbors, (row-1), col+0);
			add_neighbors(&neighbors, (row-1), col+1);
			add_neighbors(&neighbors, (row+0), col+1);
			add_neighbors(&neighbors, (row+1), col+1);

int neighbor_cnt = bin_cnt[(row-1)*bin_dim + col-1] +
bin_cnt[(row+0)*bin_dim + col-1] +
bin_cnt[(row+1)*bin_dim + col-1] +
bin_cnt[(row+1)*bin_dim + col+0] +
bin_cnt[(row-1)*bin_dim + col+0] +
bin_cnt[(row-1)*bin_dim + col+1] +
bin_cnt[(row+0)*bin_dim + col+1] +
bin_cnt[(row+1)*bin_dim + col+1];

        for (int j = 0; j < neighbor_cnt; ++j) {
            apply_force(parts[i], neighbors[j]);
        }
    }


    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}
