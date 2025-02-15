#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>

#include "../include/types.h"


void parseParticuleFile(Particles& particules, const std::string &filename, int N) {
  std::string line;
  std::ifstream particules_file(filename);
  u32 i = 0;

  if (!particules_file) {
    throw std::runtime_error("[ERROR] While opening file  " + filename);
  }

  // Ignore first line
  std::getline(particules_file, line);

  if (particules_file.is_open()) {
    while (std::getline(particules_file, line) && i < N) {
      std::istringstream lineStream(line);
      int ignore;

      // Default value of a particules
      f64 x, y, z;
      f64 fx = 0.0f, fy = 0.0f, fz = 0.0f;
      f64 mx = 0.0f, my = 0.0f, mz = 0.0f;

      lineStream >> ignore >> x >> y >> z;

      particules.x[i] = x;
      particules.y[i] = y;
      particules.z[i] = z;

      particules.fx[i] = fx;
      particules.fy[i] = fy;
      particules.fz[i] = fz;

      particules.mx[i] = mx;
      particules.my[i] = my;
      particules.mz[i] = mz;

      i++;
    }
    particules_file.close();
  }

  else
    std::cout << "Unable to open file";
}

void parse_command_line_arguments(int argc, char *argv[], SimulationParameters &params) {
  const struct option long_options[] = {
      {"steps", required_argument, nullptr, 'n'},
      {"timestep", required_argument, nullptr, 'd'},
      {"dt", required_argument, nullptr, 'd'},
      {"temperature", required_argument, nullptr, 'T'},
      {"input", required_argument, nullptr, 'i'},
      {"output", required_argument, nullptr, 'o'},
      {"cutoff", required_argument, nullptr, 'r'},
      {"rc", required_argument, nullptr, 'r'},
      {"boxsize", required_argument, nullptr, 'L'},
      {"mass", required_argument, nullptr, 'm'},
      {"relaxationtime", required_argument, nullptr, 't'},
      {"tau", required_argument, nullptr, 't'},
      {"particles", required_argument, nullptr, 'N'},
      {"seed", required_argument, nullptr, 's'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0} // Fin des options
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "n:d:T:i:o:r:L:m:t:N:s:h", long_options, nullptr)) != -1) {
    switch (opt) {
    case 'n':
      params.steps = std::stoi(optarg);
      break;
    case 'd':
      params.dt = std::stod(optarg);
      break;
    case 'T':
      params.temperature = std::stod(optarg);
      break;
    case 'i':
      params.input_file = optarg;
      break;
    case 'o':
      params.output_file = optarg;
      break;
    case 'r':
      params.cutoff = std::stod(optarg);
      break;
    case 'L': {
      params.box_size = std::stod(optarg);
      const float L = params.box_size;
      const std::vector<Vec3> periodic_images = {
          {0, 0, 0},  {-L, -L, 0}, {-L, -L, L}, {-L, 0, -L},  {-L, 0, 0},
          {-L, 0, L}, {-L, L, -L}, {-L, L, 0},  {-L, L, L},   {0, -L, -L},
          {0, -L, 0}, {0, -L, L},  {0, 0, -L},  {-L, -L, -L}, {0, 0, L},
          {0, L, -L}, {0, L, 0},   {0, L, L},   {L, -L, -L},  {L, -L, 0},
          {L, -L, L}, {L, 0, -L},  {L, 0, 0},   {L, 0, L},    {L, L, -L},
          {L, L, 0},  {L, L, L}};
    } break;
    case 'm':
      params.mass = std::stod(optarg);
      break;
    case 't':
      params.relaxation_time = std::stod(optarg);
      break;
    case 'N': {
      params.num_particles = std::stoi(optarg);

      if(params.num_particles > 1000) {
        std::cout << " The number of particles can't be greater than 1000 !" << std::endl;
        std::cout << "  -> You executed " << argv[0] << " -N " << params.num_particles << std::endl;
        exit(0);
      }
    } break;
    case 's':
      params.seed = std::stoi(optarg);
      break;
    case 'h':
      std::cout
          << "Usage: " << argv[0] << " [options]\n"
          << "Options:\n"
          << "  -n, --steps <value>        Number of time steps (default: "
          << params.steps << ")\n"
          << "  -dt, --timestep <value>    Time step size (default: "
          << params.dt << ")\n"
          << "  -T, --temperature <value>  Target temperature (default: "
          << params.temperature << ")\n"
          << "  -i, --input <file>         Input file with initial particle "
             "positions (default: "
          << params.input_file << ")\n"
          << "  -o, --output <file>        Output file for results (default: "
          << params.output_file << ")\n"
          << "                             Example for creating "
             "../out/ex_output.csv and ./out/ex_output.pdb:\n"
          << "                             " << argv[0] << "-o ex_output" << "\n"
          << "  -rc, --cutoff <value>      Cutoff distance for Lennard-Jones "
             "potential (default:"
          << params.cutoff << ")\n"
          << "  -L, --boxsize <value>      Simulation box size (default: "
          << params.box_size << ")\n"
          << "  -m, --mass <value>         Particle mass (default: "
          << params.mass << ")\n"
          << "  -tau, --relaxationtime <value>  Thermostat relaxation time "
             "(default: "
          << params.relaxation_time << ")\n"
          << "  -N, --particles <value>    Number of particles (default: "
          << params.num_particles << ")\n"
          << "  -s, --seed <value>         Random seed (default: "
          << params.seed << ")\n"
          << "  -h, --help                 Display this help message\n";
      exit(0);
    default:
      std::cerr << "Invalid option. Use -h or --help for usage information.\n";
      exit(1);
    }
  }
}

/**
 * @brief Creates the CSV file with the header if it does not already exist.
 */
void create_csv_file(const std::string filename) {
  // Ensure the directory exists
  std::filesystem::create_directories("../out");

  // Check if the file already exists
  if (std::filesystem::exists(filename)) {
    std::cout << "[INFO] The file \"" << filename << "\" already exists and will be erased." << std::endl;
  }

  std::ofstream file(filename);
  file << "time,temperature,kinetic_energy\n"; // Write the header
  file.close();
}


/**
 * @brief Creates the PDB file with the header if it does not already exist.
 */
void create_pdb_file(const std::string filename) {
  // Ensure the directory exists
  std::filesystem::create_directories("../out");

  // Check if the file already exists
  if (std::filesystem::exists(filename)) {
    std::cout << "[INFO] The file \"" << filename << "\" already exists and will be erased." << std::endl;
  }

  std::ofstream file(filename);
  file.close();
}

/**
 * @brief Appends a new step to a PDB file with carbon atom positions.
 *
 * This function writes a new MODEL section to the PDB file, including:
 * - The CRYST1 record with periodic box dimensions.
 * - The MODEL number corresponding to the iteration.
 * - ATOM records with carbon atom coordinates.
 * - The termination records (TER and ENDMDL).
 *
 * @param filename Path to the PDB file.
 * @param particles Vector containing carbon atom positions.
 * @param iteration Current simulation step (MODEL number in the PDB file).
 * @param box_size Size of the periodic simulation box.
 */
void append_pdb(const std::string &filename,
                const Particles &particles, int iteration,
                f64 box_size, u32 N) {
    // Open file in append mode
    std::ofstream file(filename, std::ios::app);
    if (!file) {
        throw std::runtime_error("[ERROR] Unable to open file: " + filename);
    }

    // 🔹 Write the periodic box information (CRYST1) and model iteration
    file << "CRYST1 "
         << std::setw(9) << std::fixed << std::setprecision(3) << box_size
         << std::setw(9) << box_size
         << std::setw(9) << box_size
         << "  90.00  90.00  90.00 P  1\n";

    file << "MODEL " << std::setw(8) << iteration << "\n";

    // 🔹 Write carbon atom data (ATOM records)
    for (size_t i = 0; i < N; ++i) {
        file << "ATOM  "
             << std::setw(5) << (i + 1)  // Atom serial number (1-based index)
             << "  C               "     // Atom name (C for carbon, no residue)
             << std::setw(8) << std::fixed << std::setprecision(3) << particles.x[i]
             << std::setw(8) << particles.y[i]
             << std::setw(8) << particles.z[i]
             << "  1.00  0.00          C  \n";  // Occupancy, temperature factor, element
    }

    // 🔹 Add termination markers
    file << "TER\n";
    file << "ENDMDL\n";

    file.close();
}

/**
 * @brief Appends a line to a CSV file storing temperature and kinetic energy
 * over time.
 *
 * If the file does not exist, it creates it and adds the header.
 *
 * @param time Current simulation time
 * @param temperature Current temperature
 * @param kinetic_energy Current kinetic energy
 */
void append_csv(const std::string filename, f64 time, f64 temperature, f64 kinetic_energy) {
  // Ensure the directory exists
  std::filesystem::create_directories("../out");

  // Open the file in append mode
  std::ofstream file(filename, std::ios::app);

  // Check if the file is empty (first write), and add headers if needed
  if (file.tellp() == 0) {
    file << "time,temperature,kinetic_energy\n";
  }

  // Append the new line with simulation data
  file << time << "," << temperature << "," << kinetic_energy << "\n";

  file.close();
}

/**
 * @brief Displays a summary of the Molecular Dynamics simulation results.
 *
 * @param boundary_conditions Indicates whether periodic boundary conditions
 * (PBC) are used.
 * @param total_simulation_time Total elapsed simulation time (in picoseconds).
 * @param computation_time Total computation time (in seconds).
 * @param num_threads Number of OpenMP threads used.
 * @param final_temperature Final computed temperature (Kelvin).
 * @param kinetic_energy Final kinetic energy of the system.
 * @param potential_energy Final potential energy.
 * @param total_energy Total system energy.
 * @param sum_forces Magnitude of the sum of all forces (should be close to 0).
 * @param total_momentum Magnitude of total system momentum (should be close to
 * 0).
 * @param params Give some informations like the number of particules, steps, box 
 * size and dt.
 */
void display_simulation_summary(bool boundary_conditions,
                                f64 total_simulation_time, int num_threads,
                                f64 final_temperature, f64 kinetic_energy,
                                f64 potential_energy, f64 total_energy,
                                f64 totalFx, f64 totalFy, f64 totalFz,
                                f64 total_momentum,
                                SimulationParameters params) {
  std::cout << std::endl << "===========================================" << std::endl;
  std::cout << " Molecular Dynamics Simulation Summary " << std::endl;
  std::cout << "===========================================" << std::endl;

  // 🔹 System Properties
  std::cout << " 🔹 Number of particles        : " << params.num_particles << std::endl;
  std::cout << " 🔹 Steps                      : " << params.steps << std::endl;
  std::cout << " 🔹 Simulation box size (L)    : " << params.box_size << "^3" << std::endl;
  std::cout << " 🔹 Boundary conditions        : "
            << (boundary_conditions ? "Periodic (27)" : "None") << std::endl;
  std::cout << " 🔹 Time step (dt)             : " << params.dt << " fs" << std::endl;
  std::cout << " 🔹 Total simulation time      : " << total_simulation_time
            << " sec" << std::endl;
  std::cout << " 🔹 Total simulation time/step : "
            << total_simulation_time / params.steps << " sec" << std::endl;
  std::cout << " 🔹 Number of threads (OMP)    : " << num_threads << std::endl;  

  // Thermodynamic Properties
  std::cout << std::endl << " \e[1mThermodynamic Properties:\e[0m" << std::endl;
  std::cout << " 🔹 Final temperature (Tfinal) : " << final_temperature << " K"
            << std::setw(8) << YELLOW
            << (final_temperature > 1000 ? "(!! So HOT)" : "")
            << RESET << std::endl;
  std::cout << " 🔹 Kinetic energy (Ek)        : " << kinetic_energy << std::endl;
  std::cout << " 🔹 Potential energy (Ep)      : " << potential_energy << std::endl;
  std::cout << " 🔹 Total energy (Etot)        : " << total_energy << std::setw(10) << YELLOW
            << (std::abs(total_energy) < 1e-3 ? "Stable" : "(!! Possible Drift)")
            << RESET << std::endl;

  // Forces & Momentum
  std::cout << std::endl << " \e[1mForces & Momentum:\e[0m" << std::endl;
  std::cout << " 🔹 Sum of forces (∑Fx)        : " << totalFx << std::setw(11) << YELLOW
            << (totalFx < 1e-3 ? "" : "(!! Check System)") << RESET
            << std::endl;
  std::cout << " 🔹 Sum of forces (∑Fy)        : " << totalFy << std::setw(10) << YELLOW
            << (totalFy < 1e-3 ? "" : "(!! Check System)") << RESET
            << std::endl;
  std::cout << " 🔹 Sum of forces (∑Fz)        : " << totalFz << std::setw(10) << YELLOW
            << (totalFz < 1e-3 ? "" : "(!! Check System)") << RESET
            << std::endl;
  std::cout << " 🔹 Total momentum (P)         : " << total_momentum << std::setw(10) << YELLOW
            << (total_momentum < 1e-3 ? "" : "(!! Possible Drift)") << RESET 
            << std::endl;

  std::cout << std::endl << BOLDGREEN << " Simulation Completed Successfully" << RESET << std::endl << std::endl;
}