#pragma once

#include <string>
#include <vector>
#include "types.h"

void parseParticuleFile(Particles& particules, const std::string &filename,
                        int N_particules_total);

void create_csv_file(const std::string filename);
void create_pdb_file(const std::string filename);

void append_csv(const std::string filename, double time, double temperature,
                   double kinetic_energy);
void append_pdb(const std::string &filename,
                const Particles &particles, int iteration,
                double box_size);

std::string getInputFilePath(const std::string &filename);

void display_simulation_summary(int num_particles, f64 box_size,
                                bool boundary_conditions, f64 dt,
                                f64 total_simulation_time, int num_threads,
                                f64 final_temperature, f64 kinetic_energy,
                                f64 potential_energy, f64 total_energy,
                                f64 totalFx, f64 totalFy, f64 totalFz, f64 total_momentum);