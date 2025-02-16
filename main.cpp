#include <chrono>
#include <ctime>

#include "include/types.h"
#include "include/simulation.h"
#include "include/utils.h"


int main(int argc, char* argv[]) {
  SimulationParameters params;

  parse_command_line_arguments(argc, argv, params);

  // Start - stop timer and execute the molecular dynamic simulation
  auto start = std::chrono::system_clock::now();
  Particles p = molecular_simulation(params);
  auto end = std::chrono::system_clock::now();

  // Variables to print
  f64 totalFx, totalFy, totalFz;
  f64 kinetic_energy, final_temperature;
  f64 U_total, total_energy;
  f64 total_momentum;
  std::chrono::duration<double> elapsed_seconds = end - start;
  bool boundary_conditions = params.N_sym == 1 ? false : true;

  // Compute variables to print
  compute_sum_forces(p, totalFx, totalFy, totalFz, params.num_particles);
  compute_kinetic_energy(p, kinetic_energy, params.num_particles, params.mass);
  compute_kinetic_temperature_energy(kinetic_energy, params.num_particles, final_temperature);
  compute_potential_energy(p, U_total, params.periodic_images, params.num_particles, params.cutoff, params.N_sym);
  compute_total_energy(U_total, kinetic_energy, total_energy);
  compute_total_momentum(p, total_momentum, params.num_particles, params.mass);

  // Display variables
  display_simulation_summary(boundary_conditions, elapsed_seconds.count(), 1,
                             final_temperature, kinetic_energy, U_total,
                             total_energy, totalFx, totalFy, totalFz,
                             total_momentum, params);

  return 0;
}
