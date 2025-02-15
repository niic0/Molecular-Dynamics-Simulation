#pragma once

#include "types.h"
#include <vector>

// 🔹 Physical parameters
constexpr f64 epsilon = 0.2;               // Depth of the potential well (in LJ units)
constexpr f64 r_star = 3.0;                // Equilibrium distance of the Lennard-Jones potential
constexpr f64 r_star_squared = r_star * r_star; // Square of r_star to avoid redundant calculations

// 🔹 Box and periodic boundary conditions parameters
constexpr f64 L = 42.0;                    // Box dimension (must be > R_cut)
constexpr f64 R_cut = 10.0;                // Cutoff radius (must be < L)
constexpr f64 R_cut_squared = R_cut * R_cut; // Cutoff radius squared
constexpr u32 N_sym = 27;                  // Number of vector translation

// 🔹 Simulation parameters
constexpr u32 N_particules_total = 1000;   // Total number of particles
constexpr u32 N_particules_local = 100;    // Number of local particles for a subset

// 🔹 Conversion factors and physical constants
constexpr f64 force_conversion = 0.0001 * 4.186; // Force unit conversion factor
constexpr f64 dt = 1.0;                   // Time step (1 femtosecond)
constexpr f64 mass = 18.0;                // Mass of a particle (e.g., water in atomic units)
constexpr f64 mass_inv = 1.0 / mass;      // Inverse mass for optimized calculations
constexpr f64 r_const = 0.00199;          // Gas constant in simulation units

// 🔹 Thermodynamic parameters
constexpr f64 N_dl = 3.0 * N_particules_total - 3.0; // Degrees of freedom of the system
constexpr f64 kinetic_temperature_factor = 1.0 / (N_dl * r_const); // Kinetic temperature factor
constexpr f64 kinetic_energy_factor = 1.0 / (2.0 * force_conversion); // Kinetic energy factor
constexpr f64 T0 = 300.0;                // Initial temperature (Kelvin)

constexpr f64 target_kinetic_energy = N_dl * r_const * T0; // Target kinetic energy for T0

// 🔹 Thermostat parameters
constexpr f64 berendsen_gamma = 0.01;    // Damping coefficient of the Berendsen thermostat

// Structure for periodic images
struct Vec3 {
  float x, y, z;
};

// Translation vectors
const std::vector<Vec3> periodic_images = {
  {0, 0, 0}, {-L, -L, 0}, {-L, -L, L},
  {-L,  0, -L}, {-L,  0, 0}, {-L,  0, L},
  {-L,  L, -L}, {-L,  L, 0}, {-L,  L, L},
  {0, -L, -L},  {0, -L, 0},  {0, -L, L},
  {0,  0, -L},  {-L, -L, -L},  {0,  0, L},
  {0,  L, -L},  {0,  L, 0},  {0,  L, L},
  {L, -L, -L},  {L, -L, 0},  {L, -L, L},
  {L,  0, -L},  {L,  0, 0},  {L,  0, L},
  {L,  L, -L},  {L,  L, 0},  {L,  L, L}
};