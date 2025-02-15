#pragma once

#include "types.h"
#include <vector>
#include <string>

// Floats
using f32 = float;
using f64 = double;
using f128 = long double;

// Unsigned int
using u16 = unsigned int;
using u32 = unsigned long int;
using u64 = unsigned long long int;

// Structure for periodic images
struct Vec3 {
  f64 x, y, z;
};

// Particules representation, structure of array
struct Particles {
  std::vector<f64> x, y , z; // Position
  std::vector<f64> fx, fy, fz; // Forces
  std::vector<f64> mx, my, mz; // momentum

  Particles(u32 particules_number) {
    x.reserve(particules_number);
    y.reserve(particules_number);
    z.reserve(particules_number);

    fx.reserve(particules_number);
    fy.reserve(particules_number);
    fz.reserve(particules_number);

    mx.reserve(particules_number);
    my.reserve(particules_number);
    mz.reserve(particules_number);
  }
};

struct SimulationParameters {
    u32 steps = 1000;                  // Nombre de pas de temps (-n)
    f64 dt = 1;                        // Pas de temps (-dt)
    f64 temperature = 300.0;           // Température cible (-T)
    std::string input_file = "../input/particule.xyz";  // Fichier d'entrée (-i)
    std::string output_file = "../out/output_data.csv";     // Fichier de sortie (-o)
    f64 cutoff = 10;                  // Distance de coupure (-rc)
    f64 box_size = 42.0;               // Taille de la boîte de simulation (-L)
    f64 mass = 18.0;                    // Masse des particules (-m)
    u32 relaxation_time = 3;           // Temps de relaxation du thermostat (-tau)
    u32 num_particles = 1000;          // Nombre de particules (-N)
    int seed = 420;                     // Graine aléatoire (-s)
    std::vector<Vec3> periodic_images = { // periodic images
        {0, 0, 0},
        {-box_size, -box_size, 0},
        {-box_size, -box_size, box_size},
        {-box_size, 0, -box_size},
        {-box_size, 0, 0},
        {-box_size, 0, box_size},
        {-box_size, box_size, -box_size},
        {-box_size, box_size, 0},
        {-box_size, box_size, box_size},
        {0, -box_size, -box_size},
        {0, -box_size, 0},
        {0, -box_size, box_size},
        {0, 0, -box_size},
        {-box_size, -box_size, -box_size},
        {0, 0, box_size},
        {0, box_size, -box_size},
        {0, box_size, 0},
        {0, box_size, box_size},
        {box_size, -box_size, -box_size},
        {box_size, -box_size, 0},
        {box_size, -box_size, box_size},
        {box_size, 0, -box_size},
        {box_size, 0, 0},
        {box_size, 0, box_size},
        {box_size, box_size, -box_size},
        {box_size, box_size, 0},
        {box_size, box_size, box_size}};
    ;
};

// Life is not in black and white
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */