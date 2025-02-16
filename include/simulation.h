#pragma once

#include <vector>
#include <string>
#include "types.h"

Particles molecular_simulation(SimulationParameters params);

void compute_velocity_verlet(Particles &p, u32 N, f64 dt, f64 mass,
                             std::vector<Vec3>& periodic_images,
                             f64 R_cut_squared, u32 N_sym);
void compute_velocity_verlet_neighbor(Particles &p, u32 N, f64 dt, f64 mass,
                                      std::vector<Vec3> &periodic_images,
                                      f64 R_cut_squared, u32 N_sym,
                                      std::vector<std::vector<u32>> neighbor);

void compute_forces_pbc_neighbor(Particles &p,
                                 std::vector<std::vector<u32>> &neighbor, u32 N,
                                 std::vector<Vec3> &periodic_images,
                                 f64 R_cut_squared, u32 N_sym);

void build_verlet_lists(Particles &p, std::vector<std::vector<u32>> &neighbor,
                        u32 N, std::vector<Vec3> &periodic_images,
                        f64 R_cut_squared, u32 N_sym, u32 n_max_neighbor);

void compute_kinetic_energy(Particles &p, f64 &kinetic_energy, u32 N, f64 mass);
void compute_kinetic_temperature_energy(f64& kinetic_energy, u32 N, f64& temperature);
void compute_sum_forces(Particles& p, f64& totalFx, f64& totalFy, f64& totalFz, u32 N);
void compute_potential_energy(Particles& p, f64& U_total, std::vector<Vec3>& periodic_images, u32 N, f64 R_cut_squared, u32 N_sym);
void compute_total_energy(f64& U_total, f64& kinetic_energy, f64& total_energy);
void compute_total_momentum(Particles& p, f64& totalMomentum, u32 N, f64 mass);

void init_momentum(Particles& p, u32 N);
void kinetic_calibration(Particles &p, f64 target_temperature, u32 N, f64 mass);

void apply_berendsen_thermostat(Particles& p, f64 T, f64 T0, u32 N);

void fix_pbc_position(Particles& p, u32 N, f64 box_size);