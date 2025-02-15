#pragma once

#include <vector>
#include <string>
#include "types.h"

Particles molecular_simulation(std::string input_particules, u32 step);
void compute_kinetic_energy(Particles& p, f64& kinetic_energy);
void compute_kinetic_temperature_energy(f64& kinetic_energy, f64& temperature);
void compute_sum_forces(Particles& p, f64& totalFx, f64& totalFy, f64& totalFz);
void compute_potential_energy(Particles& p, f64& U_total);
void compute_total_energy(f64& U_total, f64& kinetic_energy, f64& total_energy);
void compute_total_momentum(Particles& p, f64& totalMomentum);