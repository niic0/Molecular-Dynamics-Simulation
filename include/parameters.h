#pragma once

#include "types.h"
#include <vector>

constexpr f64 epsilon = 0.2;               // Depth of the potential well (in LJ units)
constexpr f64 r_star = 3.0;                // Equilibrium distance of the Lennard-Jones potential
constexpr f64 r_star_squared = r_star * r_star; // Square of r_star to avoid redundant calculations
constexpr f64 force_conversion = 0.0001 * 4.186; // Force unit conversion factor
constexpr f64 r_const = 0.00199;           // Gas constant in simulation units
constexpr f64 kinetic_energy_factor = 1.0 / (2.0 * force_conversion); // Kinetic energy factor
constexpr f64 berendsen_gamma = 0.01;     // Damping coefficient of the Berendsen thermostat