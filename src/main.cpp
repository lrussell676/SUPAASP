/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "files.h"

#include <chrono>
#include <cstdlib>

#ifdef KOKKOS_ENABLED
#include "calc_routines_kokkos.h"
#else
#include "calc_routines.h"
#endif

int main(int argc, char* argv[]) {
   auto start = std::chrono::high_resolution_clock::now(); // Start wall time
/* ------------------------------------------------------------------------------------------------
   ---- ### DEFINE SIMULATION CONDITIONS HERE ### -------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   // Timescales
   const double dt = 0.0005;           // Timestep size
   const double t = 0.0;               // Initial time
   const double t_max = 15;            // Maximum time
   const int t_pw = 100;               // Write to file every t_pw timesteps
   // Simulation Scale and Sizing
   const int N = 100;                   // Number of atoms
   const double Lx = 100.0;            // Box length in x
   const double Ly = 100.0;            // Box length in y
   const double Lz = 100.0;            // Box length in z
   int seed = 242424;                  // Random Seed for Stochastic Forces and RNGs
   
   // Atomic, Thermodynamic, and Molecular Conditions
   const double KbT = 1.0;             // Boltzmann constant
   const double gamma = 25;            // Friction coefficient
   const double length_scale = 1.0;    // Length scale for Lennard-Jones potential
   const double energy_scale = 1.0;    // Energy scale for Lennard-Jones potential
   const double rc1 = 2.5;             // Cutoff distance for Lennard-Jones potential
   const double rc2 = 2.5;             // Overlap reject thershold for random position generation
   const double T = 1.0;               // Temperature
   const double m = 1.0;               // Mass of atoms
   // Paths to write outputted data into
   std::string output_file_path_x = \
         "../src_written_data/data_pos.txt";  // Data file output for positions
   std::string output_file_path_v = \
         "../src_written_data/data_vel.txt";  // Data file output for velocities
/* ------------------------------------------------------------------------------------------------
   ---- Declaring Necessary Data from Input Conditions --------------------------------------------
   --------------------------------------------------------------------------------------------- */
   const int n_timesteps = (t_max - t) / dt;
   std::vector<double> timesteps(n_timesteps, 0.0);
   std::vector<double> iR(N, 0.0);
   const std::array<double, 3> L = {Lx, Ly, Lz};
   std::vector<std::vector<double>> R(3, std::vector<double>(N, 0.0));
   std::vector<std::vector<double>> V(3, std::vector<double>(N, 0.0));
   // Check if seed is provided as a command-line argument
   if (argc > 1) {
      seed = std::atoi(argv[1]);
      std::cout << 
      "Seed detected as command-line argument.\nUsing provided CLI seed: " << seed << std::endl;
   } else {
      std::cout << \
      "Seed not provided as a command-line argument.\nUsing default seed: " << seed << std::endl;
   }
   // Append "_RNGseed" and the seed value to file_write_paths
   output_file_path_x += "_RNGseed" + std::to_string(seed);
   output_file_path_v += "_RNGseed" + std::to_string(seed);
   // Resolve "calc_functions" class
   calc_routines* CR = NULL;
   #ifdef KOKKOS_ENABLED
   Kokkos::initialize(argc,argv);
   {
   CR = new calc_routines_kokkos;
   #else
   CR = new calc_routines;
   #endif
/* ------------------------------------------------------------------------------------------------
   ---- Calling Initialisation Functions ----------------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   std::cout << "|| --------- Simulation and Particle Initialisation ----------- ||" << std::endl;
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   CR->initialise_positions(R, N, L, rc2);
   CR->initialise_velocities(V, N, L, T, m, seed);
   CR->pbc(R, iR, L, N);
   printf("Periodic Boundary Conditions now applied to all particles.\n");
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   printf("Specified data files will be prepared, ready for new data.\n");
   prep_data_files(output_file_path_x, N, n_timesteps, L);
   prep_data_files(output_file_path_v, N, n_timesteps, L);
/* ------------------------------------------------------------------------------------------------
   ---- Performing Verlet Integration -------------------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   std::cout << "|| ------------- Beginning Verlet Integration ----------------- ||" << std::endl;
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   CR->Verlet_Integration(R, V, m, N, dt, t, t_max, n_timesteps, timesteps, \
                      output_file_path_x, output_file_path_v, t_pw,\
                      gamma, KbT, energy_scale, length_scale, rc1, L, iR, seed);
   printf("\n");
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   std::cout << "|| ----------- Successful Exit to End of Program -------------- ||" << std::endl;
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   #ifdef KOKKOS_ENABLED 
   }
   Kokkos::finalize();
   #endif
/* ------------------------------------------------------------------------------------------------
   ---- Total Wall Time ---------------------------------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> duration = end - start;
   std::cout << "Total wall time: " << duration.count() << " seconds" << std::endl;
   
   return 0;
}