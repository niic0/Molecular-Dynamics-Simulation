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
#include <iostream>

#include "../include/types.h"
#include "../include/utils.h"
#include "../include/parameters.h"

auto sign = [](auto a, auto b) { return b < 0 ? -a : a; };

/**
 * @brief Computes the forces using the Lennard-Jones potential with periodic
 * boundary conditions.
 *
 * @param p Vector of particles.
 */
void compute_forces_pbc(Particles &p) {
  for (u32 i_sym = 0; i_sym < N_sym; i_sym++) {
    for (u32 i = 0; i < N_particules_total; i++) {
      for (u32 j = i + 1; j < N_particules_total; j++) {
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
        f64 rij_inv_squared = 1 / rij_squared;

        if (rij_squared < R_cut_squared) {
          // Compute repulsion and dispersion
          f64 potential  = r_star_squared * rij_inv_squared; // (r*^2 / rij^2)
          f64 potential2 = potential * potential;            // (r*^2 / rij^2)^2
          f64 potential3 = potential2 * potential;           // (r*^2 / rij^2)^3
          f64 potential4 = potential2 * potential2;          // (r*^2 / rij^2)^4
          f64 potential7 = potential4 * potential3;          // (r*^3 / rij^4)^7

          // Compute forces
          f64 grad_inner = potential7 - potential4;
          f64 common_factor = -48 * epsilon * grad_inner;
          f64 fx = common_factor * dx;
          f64 fy = common_factor * dy;
          f64 fz = common_factor * dz;

          // Update forces on particule i and j
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

void fix_pbc_position(Particles& p) {
  f64 aint_xi_on_half_L, aint_yi_on_half_L, aint_zi_on_half_L;
  f64 half_L = L / 2.0;

  for (u32 i = 0; i < N_particules_total; i++) {
    modf(p.x[i] / half_L, &aint_xi_on_half_L);
    modf(p.y[i] / half_L, &aint_yi_on_half_L);
    modf(p.z[i] / half_L, &aint_zi_on_half_L);

    p.x[i] -= aint_xi_on_half_L * L;
    p.y[i] -= aint_yi_on_half_L * L;
    p.z[i] -= aint_zi_on_half_L * L;
  }
}

void compute_total_momentum(Particles& p, f64& totalMomentum) {
  double Px = 0.0, Py = 0.0, Pz = 0.0;

  for(u32 i = 0; i < N_particules_total; i++) {
    Px += mass * p.mx[i];
    Py += mass * p.my[i];
    Pz += mass * p.mz[i];
  }

  totalMomentum = Px + Py + Pz;
}

void compute_total_energy(f64& U_total, f64& kinetic_energy, f64& total_energy) {
  total_energy = U_total + kinetic_energy;
}

/**
   Compute kinetic energy using the formula (1/2)*force_conversion * ∑(pi^2/mi)
*/
void compute_kinetic_energy(Particles& p, f64& kinetic_energy) {
  f64 kinetic_enery_sum = 0;

  for (u32 i = 0; i < N_particules_total; i++) {
    kinetic_enery_sum += (p.mx[i] * p.mx[i] +
                          p.my[i] * p.my[i] +
                          p.mz[i] * p.mz[i]) *
                         mass_inv;
  }

  kinetic_energy = kinetic_energy_factor * kinetic_enery_sum;
}

/**
   Compute kinetic temperature energy using the formula 1/(N_dl * r_const) * kinetic_energy
*/
void compute_kinetic_temperature_energy(f64& kinetic_energy, f64& temperature) {
  temperature = kinetic_temperature_factor * kinetic_energy;
}

/**
   Recalibrate momentums to correspond the chosen temperature (300K)
*/
void kinetic_calibration(Particles& p) {
  // Compute the kinetic energy to make the momentum correspond to the temperature
  f64 kinetic_energy_init;
  compute_kinetic_energy(p, kinetic_energy_init);

  // Prevent division  by 0
  assert(kinetic_energy_init != 0 && "Division by 0");

  // Calibrate the momentum so he is coherent with TO
  f64 momentum_correction_factor = target_kinetic_energy / kinetic_energy_init;

  for(u32 i = 0; i<N_particules_total; i++) {
    p.mx[i] *= std::sqrt(momentum_correction_factor);
    p.my[i] *= std::sqrt(momentum_correction_factor);
    p.mz[i] *= std::sqrt(momentum_correction_factor);
  }
}

/**
 * @brief Computes velocity and position updates using the Velocity Verlet algorithm.
 * Forces are calculated using the Lennard-Jones potential with periodic boundary conditions.
 *
 * @param p Particles structure.
 */
void compute_velocity_verlet(Particles &p) {
  // Intermediate update of momentum (p = p - 0.5 * F * dt)
  for (u32 i = 0; i<N_particules_total; i++) {
    p.mx[i] -= 0.5 * p.fx[i] * force_conversion * dt;
    p.my[i] -= 0.5 * p.fy[i] * force_conversion * dt;
    p.mz[i] -= 0.5 * p.fz[i] * force_conversion * dt;
  }

  // Update position with the formula (x = x + (p/m) * dt)
  for(u32 i = 0; i<N_particules_total; i++) {
    p.x[i] += p.mx[i] * mass_inv * dt;
    p.y[i] += p.my[i] * mass_inv * dt;
    p.z[i] += p.mz[i] * mass_inv * dt;
  }

  // Compute forces with Lennard Jones at t + dt
  compute_forces_pbc(p);

  // Final update of momentum (p = p - 0.5 * F * dt)
  for(u32 i = 0; i<N_particules_total; i++) {
    p.mx[i] -= 0.5 * p.fx[i] * force_conversion * dt;
    p.my[i] -= 0.5 * p.fy[i] * force_conversion * dt;
    p.mz[i] -= 0.5 * p.fz[i] * force_conversion * dt;
  }
}

/**
   Init momentum with random values
*/
void init_momentum(Particles& p) {
  srand(420);

  // Ensure a balance distribution of momentums
  for (u32 i = 0; i < N_particules_total; i++) {
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
 * @brief Applies the Berendsen thermostat to scale momenta according to the target temperature.
 *
 * This function scales the momenta using the Berendsen factor:
 *      m_i = m_i + gamma * ((T0/T) -1) * m_i
 *
 * @param p Vector of particles.
 * @param T Current kinetic temperature.
 */
void apply_berendsen_thermostat(Particles& p, f64 T) {
  assert (T != 0 && "Division by 0");
  f64 thermostat_factor = berendsen_gamma * ((T0 / T) - 1.0);

  for (u32 i = 0; i<N_particules_total; i++) {
    p.mx[i] += thermostat_factor * p.mx[i];
    p.my[i] += thermostat_factor * p.my[i];
    p.mz[i] += thermostat_factor * p.mz[i];
  }
}

/**
 * @brief Computes the total potential energy of the system using the Lennard-Jones potential.
 *
 * @param particles Particles structure containing their positions.
 * @return Total potential energy of the system.
 */
void compute_potential_energy(Particles& p, f64& U_total) {
  f64 U_lj;

  U_total = 0, U_lj = 0;

  for (u32 i_sym = 0; i_sym < N_sym; i_sym++) {
    for (u32 i = 0; i < N_particules_total; i++) {
      for (u32 j = i + 1; j < N_particules_total; j++) {
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
        f64 rij_inv_squared = 1 / rij_squared;

        if (rij_squared < R_cut_squared) {
          // Compute repulsion and dispersion
          f64 potential = r_star_squared * rij_inv_squared; // (r*^2 / rij^2)
          f64 potential2 = potential * potential;           // (r*^2 / rij^2)^2
          f64 potential3 = potential2 * potential;          // (r*^2 / rij^2)^3
          f64 potential4 = potential2 * potential2;         // (r*^2 / rij^2)^4
          f64 potential7 = potential4 * potential3;         // (r*^3 / rij^4)^7

          // Compute forces
          f64 grad_inner = potential7 - potential4;
          f64 U_lj = -48 * epsilon * grad_inner;
          U_total += U_lj;
        }
      }
    }
  }
}


/**
 * @brief Compute the total energy of the system for each directions.
 * 
 * @param totalFx Total energy for the X axe
 * @param totalFy Total energy for the X axe
 * @param totalFz Total energy for the X axe
*/
void compute_sum_forces(Particles& p, f64& totalFx, f64& totalFy, f64& totalFz) {
  totalFx = 0.0f, totalFy = 0.0f, totalFz = 0.0f;

  for (u32 i = 0; i < N_particules_total; i++) {
    totalFx += p.fx[i];
    totalFy += p.fy[i];
    totalFz += p.fz[i];
  }
}

/**
   Whole simulation
*/
Particles molecular_simulation(std::string input_particules, u32 step) {
  // Input and Output filename
  std::string csv_output_filename = "../out/particules_data.csv";
  std::string pdb_output_filename = "../out/particules_data.pdb";
  std::string particules_filename = "../input/particule.xyz";

  // Create the output csv and pdb file
  create_csv_file(csv_output_filename);
  create_pdb_file(pdb_output_filename);

  // Create particles
  Particles p(N_particules_total);
  parseParticuleFile(p, particules_filename, N_particules_total);

  // Random momentum values for each particules
  init_momentum(p);

  // Calibrate momentum so they correspond to T0
  kinetic_calibration(p);

  // Each m_step, apply Berendsen thermostat correction
  u32 m_step = 5;

  // Store kinetic energy and temperature
  f64 kinetic_energy, T;

  // Compute all forces with pbc using velocity verlet
  for (u32 i = 0; i < step; i++) {
    // Get the kinetic energy and the tempeerature for csv output
    compute_kinetic_energy(p, kinetic_energy);
    compute_kinetic_temperature_energy(kinetic_energy, T);

    // Adding line to the csv output
    append_csv(csv_output_filename, i, T, kinetic_energy);

    // Add particules position to the pdb file
    append_pdb(pdb_output_filename, p, i, L);

    // One step of the computation
    compute_velocity_verlet(p);

    // Kinetic correction using Berendsen thermostat
    if (i % m_step == 0) {
      apply_berendsen_thermostat(p, T);
      fix_pbc_position(p);
      //kinetic_calibration(p); // bouuuuuh cheating
    }
  }

  // Add last particules position to the pdb file
  append_pdb(pdb_output_filename, p, step, L);

  return p;
}
