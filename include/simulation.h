#pragma once

#include <vector>
#include <string>
#include "types.h"

Particles molecular_simulation(SimulationParameters params);
void compute_kinetic_energy(Particles& p, f64& kinetic_energy, u32 N, f64 mass);
void compute_kinetic_temperature_energy(f64& kinetic_energy, u32 N, f64& temperature);
void compute_sum_forces(Particles& p, f64& totalFx, f64& totalFy, f64& totalFz, u32 N);
void compute_potential_energy(Particles& p, f64& U_total, std::vector<Vec3> periodic_images, u32 N, f64 R_cut_squared);
void compute_total_energy(f64& U_total, f64& kinetic_energy, f64& total_energy);
void compute_total_momentum(Particles& p, f64& totalMomentum, u32 N, f64 mass);