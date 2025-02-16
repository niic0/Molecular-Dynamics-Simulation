#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <omp.h>
#include <omp.h>
#include <vector>
#include <cassert>
#include <cmath>

#include "../include/types.h"
#include "../include/utils.h"
#include "../include/parameters.h"
#include "../include/simulation.h"

auto sign = [](auto a, auto b) { return b < 0 ? -a : a; };

/**
 * @brief Runs a molecular dynamics simulation using the provided parameters.
 *
 * This function performs a molecular dynamics simulation over a specified
 * number of time steps. It initializes the particles, computes their
 * interactions using the Velocity Verlet algorithm, and applies a Berendsen
 * thermostat to maintain the target temperature. The simulation results
 * (particle positions, temperature, and kinetic energy) are saved to CSV and
 * PDB files for analysis.
 *
 * @param params SimulationParameters object containing all simulation
 * parameters (e.g., number of particles, time step, target temperature, box
 * size, etc.).
 * @return Particles object containing the final state of the particles after
 * the simulation.
 *
 * @note The simulation uses periodic boundary conditions (PBC) and a cutoff
 * distance for the Lennard-Jones potential.
 * @note The Berendsen thermostat is applied periodically to adjust the system's
 * temperature.
 */
Particles molecular_simulation(SimulationParameters params) {
  // Construct file paths for output files
  std::string csv_output_filename = "../out/" + params.output_file + ".csv";
  std::string pdb_output_filename = "../out/" + params.output_file + ".pdb";
  std::string particules_filename = params.input_file;

  // Extract simulation parameters for convenience
  u32 N = params.num_particles;
  f64 mass = params.mass;
  u32 steps = params.steps;
  f64 target_temperature = params.temperature;
  f64 box_size = params.box_size;
  u32 relaxation_time = params.relaxation_time;
  f64 dt = params.dt;
  std::vector<Vec3> periodic_images = params.periodic_images;
  f64 R_cut_squared = params.cutoff * params.cutoff;
  u32 N_sym = params.N_sym;

  // Current kinetic energy and temperature variables
  f64 kinetic_energy, current_temperature;

  // Verlet list variables
  f64 density = N / (box_size * box_size * box_size);
  u32 n_max_neighbor = static_cast<u32>(density * 4 * M_PI * params.cutoff * params.cutoff * params.cutoff) * 2;
  std::vector<std::vector<u32>> neighbor(N, std::vector<u32>(2 * n_max_neighbor + 1));

  // Create the output csv and pdb file
  create_csv_file(csv_output_filename);
  create_pdb_file(pdb_output_filename);

  // Initialize particles from the input file
  Particles p(N);
  parseParticuleFile(p, particules_filename, N);

  // Initialize particle momenta with random values
  init_momentum(p, N);

  // Calibrate momenta to match the target temperature
  kinetic_calibration(p, target_temperature, N, mass);

  // Build the verlet (neighbors) list, this list is updated each
  // "relaxation_time"
  build_verlet_lists(p, neighbor, N, periodic_images, R_cut_squared, N_sym, n_max_neighbor);

  for (u32 i = 0; i < params.steps; i++) {
    // Compute kinetic energy and temperature
    compute_kinetic_energy(p, kinetic_energy, N, mass);
    compute_kinetic_temperature_energy(kinetic_energy, N, current_temperature);

    // Save current state to CSV and PDB files
    append_csv(csv_output_filename, i, current_temperature, kinetic_energy);
    append_pdb(pdb_output_filename, p, i, box_size, N);

    // Perform one step of the Velocity Verlet algorithm
    compute_velocity_verlet_neighbor(p, N, dt, mass, periodic_images, R_cut_squared,
                            N_sym, neighbor);

    // Apply the Berendsen thermostat periodically to maintain the target
    // temperature
    if (i % relaxation_time == 0) {
      apply_berendsen_thermostat(p, current_temperature, target_temperature,
                                 params.num_particles);
      build_verlet_lists(p, neighbor, N, periodic_images, R_cut_squared, N_sym,
                         n_max_neighbor);
    }

    // Ensure particles remain within the simulation box using periodic boundary
    // conditions
    fix_pbc_position(p, N, box_size);
  }

  // Compute the final state of kinetic energy and temperature
  compute_kinetic_energy(p, kinetic_energy, N, mass);
  compute_kinetic_temperature_energy(kinetic_energy, N, current_temperature);

  // Save the final state of the particles to the CSV and PDB file
  append_pdb(pdb_output_filename, p, steps, box_size, N);
  append_csv(csv_output_filename, steps, current_temperature, kinetic_energy);

  return p;
}

/**
 * @brief Computes the forces between particles using the Lennard-Jones
 * potential with periodic boundary conditions (PBC).
 *
 * This function calculates the forces between all pairs of particles in the
 * system, taking into account periodic boundary conditions. The forces are
 * computed using the Lennard-Jones potential, and interactions beyond the
 * cutoff distance (`R_cut_squared`) are ignored. The forces are updated in
 * place for each particle.
 *
 * @param p Reference to the Particles object containing the particle positions
 * and forces.
 * @param N Number of particles in the system.
 * @param periodic_images Vector of periodic image translations to apply for
 * PBC.
 * @param R_cut_squared Squared cutoff distance for the Lennard-Jones potential.
 *
 * @note This function modifies the forces (`fx`, `fy`, `fz`) of the particles
 * in place.
 * @note The Lennard-Jones potential is computed as:
 *       U(r) = -48 * epsilon * [(r_star/r)^14 - (r_star/r)^8]
 *       where sigma is the distance at which the potential is zero, and epsilon
 * is the depth of the potential well.
 * @note The cutoff distance is used to ignore interactions beyond a certain
 * range for efficiency.
 * @warning This function is not used because it doesn't use verlet list which
 * is faster.
 */
void compute_forces_pbc(Particles &p, u32 N, std::vector<Vec3> &periodic_images,
                        f64 R_cut_squared, u32 N_sym) {
  for (u32 i_sym = 0; i_sym < N_sym; i_sym++) {
    for (u32 i = 0; i < N; i++) {
      for (u32 j = i + 1; j < N; j++) {
        // Translating particule j position with periodic boundary conditions
        f64 xj_loc = p.x[j] + periodic_images[i_sym].x;
        f64 yj_loc = p.y[j] + periodic_images[i_sym].y;
        f64 zj_loc = p.z[j] + periodic_images[i_sym].z;

        // Compute squared distance between particles i and j (r_ij^2)
        f64 dx = p.x[i] - xj_loc;
        f64 dy = p.y[i] - yj_loc;
        f64 dz = p.z[i] - zj_loc;
        f64 rij_squared = dx * dx + dy * dy + dz * dz;
        f64 rij_inv_squared = 1 / rij_squared;

        if (rij_squared < R_cut_squared) {
          // Compute repulsion and dispersion
          f64 potential  = r_star_squared * rij_inv_squared; // (r*^2 / rij^2)
          f64 potential2 = potential * potential;            // (r*^2 / rij^2)^2
          f64 potential3 = potential2 * potential;           // (r*^2 / rij^2)^3
          f64 potential4 = potential2 * potential2;          // (r*^2 / rij^2)^4
          f64 potential7 = potential4 * potential3;          // (r*^2 / rij^2)^7

          // Compute forces using the gradient of the Lennard-Jones potential
          f64 grad_inner = potential7 - potential4;
          f64 common_factor = -48.0 * epsilon * grad_inner;
          f64 fx = common_factor * dx;
          f64 fy = common_factor * dy;
          f64 fz = common_factor * dz;

          // Update forces on particule i and j (Newton's third law)
          p.fx[i] += fx;
          p.fy[i] += fy;
          p.fz[i] += fz;
          p.fx[j] -= fx;
          p.fy[j] -= fy;
          p.fz[j] -= fz;
        }
      }
    }
  }
}

/**
 * @brief Performs a single time step of the Velocity Verlet integration
 * algorithm.
 *
 * This function updates the positions and momentum of all particles in the
 * system using the Velocity Verlet algorithm. The algorithm consists of three
 * main steps:
 * 1. Intermediate update of momentum using the current forces.
 * 2. Update of positions using the intermediate momentum.
 * 3. Recalculation of forces at the new positions and final update of momentum.
 *
 * @param p Reference to the Particles object containing the particle positions,
 * momentum, and forces.
 * @param N Number of particles in the system.
 * @param dt Time step for the integration.
 * @param mass Mass of the particles (assumed to be the same for all particles).
 * @param periodic_images Vector of periodic image translations to apply for
 * periodic boundary conditions (PBC).
 * @param R_cut_squared Squared cutoff distance for the Lennard-Jones potential.
 *
 * @note The forces are recalculated after updating the positions to account for
 * interactions at the new positions.
 * @note The forces are reset to zero before recalculating them to avoid
 * accumulating contributions from previous steps.
 * @warning This function is not use because it doesn't uses verlet lists so it
 * is less optimize
 */
void compute_velocity_verlet(Particles &p, u32 N, f64 dt, f64 mass,
                             std::vector<Vec3>& periodic_images,
                             f64 R_cut_squared, u32 N_sym) {
  const f64 mass_inv = 1 / mass;

  // Intermediate update of momentum (p = p - 0.5 * F * dt)
  for (u32 i = 0; i<N; i++) {
    p.mx[i] -= 0.5 * p.fx[i] * force_conversion * dt;
    p.my[i] -= 0.5 * p.fy[i] * force_conversion * dt;
    p.mz[i] -= 0.5 * p.fz[i] * force_conversion * dt;
  }

  // Update position with the formula (x = x + (p/m) * dt)
  for(u32 i = 0; i<N; i++) {
    p.x[i] += p.mx[i] * mass_inv * dt;
    p.y[i] += p.my[i] * mass_inv * dt;
    p.z[i] += p.mz[i] * mass_inv * dt;
  }

  // Reset forces to zero before computing new forces
  for (u32 i = 0; i < N; i++) {
    p.fx[i] = 0.0;
    p.fy[i] = 0.0;
    p.fz[i] = 0.0;
  }

  // Compute forces with Lennard-Jones potential at the new positions (t + dt)
  // without verlet list
  compute_forces_pbc(p, N, periodic_images, R_cut_squared, N_sym);

  // Final update of momentum (p = p - 0.5 * F * dt)
  for(u32 i = 0; i<N; i++) {
    p.mx[i] -= 0.5 * p.fx[i] * force_conversion * dt;
    p.my[i] -= 0.5 * p.fy[i] * force_conversion * dt;
    p.mz[i] -= 0.5 * p.fz[i] * force_conversion * dt;
  }
}

/**
 * @brief Builds the Verlet neighbor lists for all particles.
 *
 * This function constructs the Verlet neighbor lists for each particle, storing
 * the indices of neighboring particles and their periodic images within the
 * cutoff distance.
 *
 * @param p Reference to the Particles object containing the particle positions.
 * @param neighbor Reference to the 2D array storing the neighbor lists.
 * @param N Number of particles in the system.
 * @param periodic_images Vector of periodic image translations to apply for
 * PBC.
 * @param R_cut_squared Squared cutoff distance for the Lennard-Jones potential.
 * @param N_sym Number of periodic images to consider.
 * @param n_max_neighbor Maximum number of neighbors per particle.
 */
void build_verlet_lists(Particles &p, std::vector<std::vector<u32>> &neighbor,
                        u32 N, std::vector<Vec3> &periodic_images,
                        f64 R_cut_squared, u32 N_sym, u32 n_max_neighbor) {
  for (u32 i = 0; i < N; i++) {
    u32 ni = 0; // Number of neighbors for particle i
    for (u32 j = i+1; j < N; j++) {
      if (i == j)
        continue; // Skip self-interaction
      for (u32 i_sym = 0; i_sym < N_sym; i_sym++) {
        // Translating particle j position with periodic boundary conditions
        f64 xj_loc = p.x[j] + periodic_images[i_sym].x;
        f64 yj_loc = p.y[j] + periodic_images[i_sym].y;
        f64 zj_loc = p.z[j] + periodic_images[i_sym].z;

        // Compute squared distance between particles i and j
        f64 dx = p.x[i] - xj_loc;
        f64 dy = p.y[i] - yj_loc;
        f64 dz = p.z[i] - zj_loc;
        f64 rij_squared = dx * dx + dy * dy + dz * dz;

        // Store the neighbor index and periodic image if under cutoff
        if (rij_squared < R_cut_squared) {
          u32 mi = 2 * ni + 1;
          neighbor[i][mi] = j;         // Index of neighbor j
          neighbor[i][mi+1] = i_sym;   // Periodic image index
          ni++;
        }
      }
    }

    // Store the number of neighbors for particle i in [i,0]
    neighbor[i][0] = ni;
  }
}


/**
 * @brief Computes the forces between particles using the Lennard-Jones potential
 * with periodic boundary conditions (PBC) and Verlet neighbor lists.
 *
 * This function calculates the forces between particles using precomputed Verlet
 * neighbor lists, which store the indices of neighboring particles within the
 * cutoff distance. This reduces the computational cost compared to the naive O(N^2)
 * approach.
 *
 * @param p Reference to the Particles object containing the particle positions
 * and forces.
 * @param neighbor Reference to the 2D array storing the Verlet neighbor lists.
 * @param N Number of particles in the system.
 * @param periodic_images Vector of periodic image translations to apply for PBC.
 * @param R_cut_squared Squared cutoff distance for the Lennard-Jones potential.
 * @param N_sym Number of periodic images to consider.
 */
void compute_forces_pbc_neighbor(Particles &p, std::vector<std::vector<u32>> &neighbor, u32 N,
                        std::vector<Vec3> &periodic_images, f64 R_cut_squared, u32 N_sym) {
  // Loop over all particles
  for (u32 i = 0; i < N; i++) {
    // Get the number of neighbors for particle i
    u32 ni = neighbor[i][0];

    // Loop over all neighbors of particle i
    for (u32 k = 0; k < ni; k++) {
      // Get the index of the neighbor particle j and its periodic image
      u32 j = neighbor[i][2 * k + 1];     // Index of neighbor j
      u32 i_sym = neighbor[i][2 * k + 2]; // Periodic image index

      // Translating particle j position with periodic boundary conditions
      f64 xj_loc = p.x[j] + periodic_images[i_sym].x;
      f64 yj_loc = p.y[j] + periodic_images[i_sym].y;
      f64 zj_loc = p.z[j] + periodic_images[i_sym].z;

      // Compute squared distance between particles i and j
      f64 dx = p.x[i] - xj_loc;
      f64 dy = p.y[i] - yj_loc;
      f64 dz = p.z[i] - zj_loc;
      f64 rij_squared = dx * dx + dy * dy + dz * dz;

      // Prevent division by zero
      //assert(rij_squared != 0 && "Division by 0");
      f64 rij_inv_squared = 1.0 / rij_squared;

      // Compute repulsion and dispersion terms of the Lennard-Jones potential
      f64 potential = r_star_squared * rij_inv_squared; // (r*^2 / rij^2)
      f64 potential2 = potential * potential;           // (r*^2 / rij^2)^2
      f64 potential3 = potential2 * potential;          // (r*^2 / rij^2)^3
      f64 potential4 = potential2 * potential2;         // (r*^2 / rij^2)^4
      f64 potential7 = potential4 * potential3;         // (r*^3 / rij^4)^7

      // Compute forces using the gradient of the Lennard-Jones potential
      f64 grad_inner = potential7 - potential4;
      f64 common_factor = -48.0 * epsilon * grad_inner;
      f64 fx = common_factor * dx;
      f64 fy = common_factor * dy;
      f64 fz = common_factor * dz;

      // Update forces on particles i and j (Newton's third law)
      p.fx[i] += fx;
      p.fy[i] += fy;
      p.fz[i] += fz;
      p.fx[j] -= fx;
      p.fy[j] -= fy;
      p.fz[j] -= fz;
    }
  }
}

/**
 * @brief Performs a single time step of the Velocity Verlet integration
 * algorithm.
 *
 * This function updates the positions and momentum of all particles in the
 * system using the Velocity Verlet algorithm. The algorithm consists of three
 * main steps:
 * 1. Intermediate update of momentum using the current forces.
 * 2. Update of positions using the intermediate momentum.
 * 3. Recalculation of forces at the new positions and final update of momentum.
 *
 * @param p Reference to the Particles object containing the particle positions,
 * momentum, and forces.
 * @param N Number of particles in the system.
 * @param dt Time step for the integration.
 * @param mass Mass of the particles (assumed to be the same for all particles).
 * @param periodic_images Vector of periodic image translations to apply for
 * periodic boundary conditions (PBC).
 * @param R_cut_squared Squared cutoff distance for the Lennard-Jones potential.
 *
 * @note The forces are recalculated after updating the positions to account for
 * interactions at the new positions.
 * @note The forces are reset to zero before recalculating them to avoid
 * accumulating contributions from previous steps.
 */
void compute_velocity_verlet_neighbor(Particles &p, u32 N, f64 dt, f64 mass,
                             std::vector<Vec3>& periodic_images,
                             f64 R_cut_squared, u32 N_sym, std::vector<std::vector<u32>> neighbor) {
  const f64 mass_inv = 1 / mass;

  // Intermediate update of momentum (p = p - 0.5 * F * dt)
  for (u32 i = 0; i<N; i++) {
    p.mx[i] -= 0.5 * p.fx[i] * force_conversion * dt;
    p.my[i] -= 0.5 * p.fy[i] * force_conversion * dt;
    p.mz[i] -= 0.5 * p.fz[i] * force_conversion * dt;
  }

  // Update position with the formula (x = x + (p/m) * dt)
  for(u32 i = 0; i<N; i++) {
    p.x[i] += p.mx[i] * mass_inv * dt;
    p.y[i] += p.my[i] * mass_inv * dt;
    p.z[i] += p.mz[i] * mass_inv * dt;
  }

  // Reset forces to zero before computing new forces
  for (u32 i = 0; i < N; i++) {
    p.fx[i] = 0.0;
    p.fy[i] = 0.0;
    p.fz[i] = 0.0;
  }

  // Compute forces with Lennard-Jones potential at the new positions (t + dt)
  compute_forces_pbc_neighbor(p, neighbor, N, periodic_images, R_cut_squared, N_sym);

  // Final update of momentum (p = p - 0.5 * F * dt)
  for(u32 i = 0; i<N; i++) {
    p.mx[i] -= 0.5 * p.fx[i] * force_conversion * dt;
    p.my[i] -= 0.5 * p.fy[i] * force_conversion * dt;
    p.mz[i] -= 0.5 * p.fz[i] * force_conversion * dt;
  }
}

/**
 * @brief Computes the total potential energy of the system using the
 * Lennard-Jones potential.
 *
 * This function calculates the total potential energy of the system by summing
 * up the contributions from all pairs of particles interacting via the
 * Lennard-Jones potential. Periodic boundary conditions are applied using the
 * provided periodic images. Interactions beyond the cutoff distance
 * (`R_cut_squared`) are ignored for efficiency.
 *
 * @param p Reference to the Particles object containing the particle positions.
 * @param U_total Reference to the variable where the total potential energy
 * will be stored.
 * @param periodic_images Vector of periodic image translations to apply for
 * periodic boundary conditions (PBC).
 * @param N Number of particles in the system.
 * @param R_cut_squared Squared cutoff distance for the Lennard-Jones potential.
 *
 * @note The Lennard-Jones potential is computed as:
 *       U(r) = 4 * epsilon * [(r_star/r)^12 - (r_star/r)^6]
 *       where r_star is the distance at which the potential is zero, and
 * epsilon is the depth of the potential well.
 * @note The function ensures that the distance between particles is non-zero to
 * avoid division by zero.
 * @note The cutoff distance is used to ignore interactions beyond a certain
 * range for efficiency.
 */
void compute_potential_energy(Particles &p, f64 &U_total,
                              std::vector<Vec3>& periodic_images, u32 N,
                              f64 R_cut_squared, u32 N_sym) {
  U_total = 0.0;

  for (u32 i_sym = 0; i_sym < N_sym; i_sym++) {
    for (u32 i = 0; i < N; i++) {
      for (u32 j = i + 1; j < N; j++) {
        // Translating particule j position with periodic boundary conditions
        f64 xj_loc = p.x[j] + periodic_images[i_sym].x;
        f64 yj_loc = p.y[j] + periodic_images[i_sym].y;
        f64 zj_loc = p.z[j] + periodic_images[i_sym].z;

        // Compute r_ij^2
        f64 dx = p.x[i] - xj_loc;
        f64 dy = p.y[i] - yj_loc;
        f64 dz = p.z[i] - zj_loc;
        f64 rij_squared = dx * dx + dy * dy + dz * dz;

        // Prevent Division by 0
        assert(rij_squared != 0 && "Division by 0");
        f64 rij_inv_squared = 1.0 / rij_squared;

        if (rij_squared < R_cut_squared) {
          // Compute repulsion and dispersion
          f64 potential = r_star_squared * rij_inv_squared; // (r*^2 / rij^2)
          f64 potential2 = potential * potential;           // (r*^2 / rij^2)^2
          f64 potential3 = potential2 * potential;          // (r*^2 / rij^2)^3
          f64 potential6 = potential3 * potential3;         // (r*^2 / rij^2)^6

          // Compute forces
          f64 grad_inner = potential6 - (2 * potential3);
          f64 U_lj = 4.0 * epsilon * grad_inner;
          U_total += U_lj;
        }
      }
    }
  }
}

/**
 * @brief Computes the total force acting on the system in each direction (x, y,
 * z).
 *
 * This function calculates the sum of the forces acting on all particles in the
 * system for each spatial direction (x, y, z). The results are stored in the
 * provided reference variables `totalFx`, `totalFy`, and `totalFz`.
 *
 * @param p Reference to the Particles object containing the particle forces.
 * @param totalFx Reference to the variable where the total force in the
 * x-direction will be stored.
 * @param totalFy Reference to the variable where the total force in the
 * y-direction will be stored.
 * @param totalFz Reference to the variable where the total force in the
 * z-direction will be stored.
 * @param N Number of particles in the system.
 *
 * @note This function is useful for verifying that the net force on the system
 * is zero, which is expected in a closed system with balanced forces (e.g., due
 * to Newton's third law).
 */
void compute_sum_forces(Particles &p, f64 &totalFx, f64 &totalFy, f64 &totalFz,
                        u32 N) {
  totalFx = 0.0f, totalFy = 0.0f, totalFz = 0.0f;

  for (u32 i = 0; i < N; i++) {
    totalFx += p.fx[i];
    totalFy += p.fy[i];
    totalFz += p.fz[i];
  }
}

/**
 * @brief Calibrates the particle momentum to match a target temperature.
 *
 * This function adjusts the momenta of all particles in the system so that the
 * kinetic energy corresponds to the desired target temperature. The calibration
 * is performed by scaling the momentum of each particle by a correction factor
 * derived from the ratio of the target kinetic energy to the current kinetic
 * energy.
 *
 * @param p Reference to the Particles object containing the particle momentum.
 * @param target_temperature Desired temperature of the system.
 * @param N Number of particles in the system.
 * @param mass Mass of the particles (assumed to be the same for all particles).
 *
 * @note The target kinetic energy is computed using the equipartition theorem:
 *       E_kinetic = (N_dl * k_B * T) / 2
 *       where N_dl is the number of degrees of freedom (3N - 3 for 3D
 * translational motion), k_B is the Boltzmann constant (represented here by
 * `r_const`), and T is the target temperature.
 * @note The function assumes that the initial kinetic energy is non-zero to
 * avoid division by zero.
 */
void kinetic_calibration(Particles &p, f64 target_temperature, u32 N,
                         f64 mass) {
  // Compute the number of degrees of freedom (3N - 3 for 3D translational
  // motion)
  f64 N_dl = 3.0 * N - 3.0;
  f64 target_kinetic_energy = N_dl * r_const * target_temperature;

  // Compute the initial kinetic energy of the system
  f64 kinetic_energy_init;
  compute_kinetic_energy(p, kinetic_energy_init, N, mass);

  // Prevent division  by 0
  assert(kinetic_energy_init != 0 && "Division by 0");

  // Calibrate the momentum so he is coherent with TO
  f64 momentum_correction_factor = target_kinetic_energy / kinetic_energy_init;

  // Adjust the momenta of all particles to match the target kinetic energy
  for(u32 i = 0; i<N; i++) {
    p.mx[i] *= std::sqrt(momentum_correction_factor);
    p.my[i] *= std::sqrt(momentum_correction_factor);
    p.mz[i] *= std::sqrt(momentum_correction_factor);
  }
}

/**
 * @brief Adjusts particle positions to enforce periodic boundary conditions
 * (PBC).
 *
 * This function ensures that particles remain within the simulation box by
 * wrapping their positions using periodic boundary conditions. If a particle
 * moves outside the box, it is translated back into the box by subtracting the
 * appropriate multiple of the box size.
 *
 * @param p Reference to the Particles object containing the particle positions.
 * @param N Number of particles in the system.
 * @param box_size Size of the simulation box (assumed to be cubic).
 */
void fix_pbc_position(Particles& p, u32 N, f64 box_size) {
  f64 aint_xi_on_half_L, aint_yi_on_half_L, aint_zi_on_half_L;
  f64 half_L = box_size / 2.0;

  for (u32 i = 0; i < N; i++) {
    // Compute the integer part of the position divided by half the box size for
    // each coordinate
    modf(p.x[i] / half_L, &aint_xi_on_half_L);
    modf(p.y[i] / half_L, &aint_yi_on_half_L);
    modf(p.z[i] / half_L, &aint_zi_on_half_L);

    // Adjust the particle positions to enforce periodic boundary conditions
    p.x[i] -= aint_xi_on_half_L * box_size;
    p.y[i] -= aint_yi_on_half_L * box_size;
    p.z[i] -= aint_zi_on_half_L * box_size;
  }
}

/**
 * @brief Computes the total momentum of the system.
 *
 * This function calculates the total momentum of all particles in the system by
 * summing up the momentum components (Px, Py, Pz) for each particle. The total
 * momentum is stored in the provided reference variable `totalMomentum`.
 *
 * @param p Reference to the Particles object containing the particle momentum
 * components.
 * @param totalMomentum Reference to a variable where the total momentum will be
 * stored.
 * @param N Number of particles in the system.
 * @param mass Mass of the particles (assumed to be the same for all particles).
 *
 * @note The momentum components (mx, my, mz) in the Particles object are
 * assumed to represent the velocity components multiplied by the mass.
 */
void compute_total_momentum(Particles &p, f64 &totalMomentum, u32 N, f64 mass) {
  double Px = 0.0, Py = 0.0, Pz = 0.0;

  for (u32 i = 0; i < N; i++) {
    Px += mass * p.mx[i];
    Py += mass * p.my[i];
    Pz += mass * p.mz[i];
  }

  // Compute the total momentum as the norm of thesum of the components
  totalMomentum = std::sqrt(Px * Px + Py * Py + Pz * Pz);
}

/**
 * @brief Computes the total kinetic energy of the system using using
 * (1/2) * force_conversion * ∑(pi^2/mi)
 *
 * This function calculates the total kinetic energy of all particles in the
 * system by summing up the kinetic energy contributions from each particle. The
 * kinetic energy for each particle is computed using its momentum components
 * (mx, my, mz) and the particle mass. The result is stored in the provided
 * reference variable `kinetic_energy`.
 *
 * @param p Reference to the Particles object containing the particle momentum
 * components.
 * @param kinetic_energy Reference to the variable where the total kinetic
 * energy will be stored.
 * @param N Number of particles in the system.
 * @param mass Mass of the particles (assumed to be the same for all particles).
 *
 * @note The momentum components (mx, my, mz) in the Particles object are
 * assumed to represent the velocity components multiplied by the mass.
 * @note The kinetic energy is scaled by a global factor
 * `kinetic_energy_factor`.
 */
void compute_kinetic_energy(Particles& p, f64& kinetic_energy, u32 N, f64 mass) {
  f64 kinetic_enery_sum = 0;
  f64 mass_inv = 1.0 / mass;

  for (u32 i = 0; i < N; i++) {
    kinetic_enery_sum +=
        (p.mx[i] * p.mx[i] + p.my[i] * p.my[i] + p.mz[i] * p.mz[i]) * mass_inv;
  }

  // Scale the total kinetic energy by the global factor
  kinetic_energy = kinetic_energy_factor * kinetic_enery_sum;
}

/**
 * @brief Computes the total energy of the system.
 *
 * This function calculates the total energy of the system by summing the
 * potential energy (U_total) and the kinetic energy (kinetic_energy). The
 * result is stored in the provided reference variable `total_energy`.
 *
 * @param U_total Reference to the variable containing the total potential
 * energy of the system.
 * @param kinetic_energy Reference to the variable containing the total kinetic
 * energy of the system.
 * @param total_energy Reference to the variable where the total energy will be
 * stored.
 *
 * @note This function assumes that U_total and kinetic_energy have already been
 * computed and provided as inputs.
 * @warning The total energy can't be used as stability criteria as we use a
 * Berendsen thermostat ! The Berendsen thermostat modify the position to
 * correspond to a certain temperature so the total energy is not conserved.
 */
void compute_total_energy(f64& U_total, f64& kinetic_energy, f64& total_energy) {
  total_energy = U_total + kinetic_energy;
}

/**
 * @brief Computes the kinetic temperature of the system from the total kinetic
 * energy.
 *
 * This function calculates the kinetic temperature of the system using the
 * total kinetic energy. The kinetic temperature is derived from the
 * equipartition theorem, which relates the kinetic energy to the temperature of
 * the system. The result is stored in the provided reference variable
 * `temperature`.
 *
 * @param kinetic_energy Reference to the variable containing the total kinetic
 * energy of the system.
 * @param N Number of particles in the system.
 * @param temperature Reference to the variable where the computed temperature
 * will be stored.
 *
 * @note The kinetic temperature is computed using the formula:
 *       T = (2 * kinetic_energy) / (N_dl * k_B)
 *       where N_dl is the number of degrees of freedom (3N - 3 for a system
 * with 3D translational motion), and k_B is the Boltzmann constant (represented
 * here by `r_const`).
 */
void compute_kinetic_temperature_energy(f64 &kinetic_energy, u32 N,
                                        f64 &temperature) {
  // Compute the number of degrees of freedom (3N - 3 for 3D translational
  // motion)
  const f64 N_dl = 3.0 * N - 3.0;

  // Compute the factor to convert kinetic energy to temperature
  const f64 kinetic_temperature_factor = 1.0 / (N_dl * r_const);

  // Compute the kinetic temperature
  temperature = kinetic_temperature_factor * kinetic_energy;
}

/**
 * @brief Initializes the momentum of particles with random values.
 *
 * This function initializes the momentum of all particles in the system with
 * random values. The momentum are generated using a uniform random distribution
 * and are balanced to ensure that the system has no net momentum.
 *
 * @param p Reference to the Particles object where the momenta will be stored.
 * @param N Number of particles in the system.
 * @param seed The seed of the random generator.
 *
 * @note The momenta are generated such that they are symmetrically distributed
 * around zero, ensuring that the total momentum of the system is approximately
 * zero.
 * @note The random number generator is seeded with a value for
 * reproducibility or diversity of generation.
 */
void init_momentum(Particles& p, u32 N) {
  srand(420);

  // Ensure a balance distribution of momentums
  for (u32 i = 0; i < N; i++) {
    // Generate a random number for each momentum direction
    f64 cx = (f64)rand() / (f64)RAND_MAX; // Random between 0 and 1
    f64 cy = (f64)rand() / (f64)RAND_MAX;
    f64 cz = (f64)rand() / (f64)RAND_MAX;
    f64 sx = (f64)rand() / (f64)RAND_MAX;
    f64 sy = (f64)rand() / (f64)RAND_MAX;
    f64 sz = (f64)rand() / (f64)RAND_MAX;

    p.mx[i] = sign(1.0, 0.5 - sx) * cx;
    p.my[i] = sign(1.0, 0.5 - sy) * cy;
    p.mz[i] = sign(1.0, 0.5 - sz) * cz;
  }
}

/**
 * @brief Applies the Berendsen thermostat to scale particle momentum.
 *
 * This function adjusts the momenta of all particles in the system using the
 * Berendsen thermostat to bring the system's temperature closer to the target
 * temperature (T0). The scaling factor is computed based on the current
 * temperature (T) and the target temperature (T0).
 *
 * @param p Reference to the Particles object containing the particle momenta.
 * @param T Current temperature of the system.
 * @param T0 Target temperature of the system.
 * @param N Number of particles in the system.
 *
 * @note The Berendsen thermostat scales the momenta by a factor that depends on
 * the difference between the current temperature (T) and the target temperature
 * (T0), as well as the relaxation parameter `berendsen_gamma`.
 * @note The function asserts that the current temperature (T) is not zero to
 * avoid division by zero.
 */
void apply_berendsen_thermostat(Particles& p, f64 T, f64 T0, u32 N) {
  assert (T != 0 && "Division by 0");
  f64 thermostat_factor = berendsen_gamma * ((T0 / T) - 1.0);

  for (u32 i = 0; i<N; i++) {
    p.mx[i] += thermostat_factor * p.mx[i];
    p.my[i] += thermostat_factor * p.my[i];
    p.mz[i] += thermostat_factor * p.mz[i];
  }
}