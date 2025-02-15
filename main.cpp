#include <chrono>
#include <ctime>

#include "include/parameters.h"
#include "include/types.h"
#include "include/simulation.h"
#include "include/utils.h"


int main() {
  std::string input_particule_file_path = "../input/particule.xyz";
  u32 step = 100;

  // Start - stop timer and execute the molecular dynamic simulation
  auto start = std::chrono::system_clock::now();
  Particles p = molecular_simulation(input_particule_file_path, step);
  auto end = std::chrono::system_clock::now();

  // Variables to display
  f64 totalFx, totalFy, totalFz;
  f64 kinetic_energy, final_temperature;
  f64 U_total, total_energy;
  f64 total_momentum;
  std::chrono::duration<double> elapsed_seconds = end - start;

  // Compute variable to display
  compute_sum_forces(p, totalFx, totalFy, totalFz);
  compute_kinetic_energy(p, kinetic_energy);
  compute_kinetic_temperature_energy(kinetic_energy, final_temperature);
  compute_potential_energy(p, U_total);
  compute_total_energy(U_total, kinetic_energy, total_energy);
  compute_sum_forces(p, totalFx, totalFy, totalFz);
  compute_total_momentum(p, total_momentum);

  // Display variables
  display_simulation_summary(N_particules_total, L, true, 1,
                                  elapsed_seconds.count(), 1, final_temperature,
                                  kinetic_energy, U_total, total_energy,
                                  totalFx, totalFy, totalFz, total_momentum);

  return 0;
}
