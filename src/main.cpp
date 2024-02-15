/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines.h"
#include "files.h"

#include <chrono>

int main() {
   auto start = std::chrono::high_resolution_clock::now(); // Start wall time
/* ------------------------------------------------------------------------------------------------
   ---- ### DEFINE SIMULATION CONDITIONS HERE ### -------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   // Timescales
   const double dt = 0.0005;           // Timestep size
   const double t = 0.0;               // Initial time
   const double t_max = 1.0;          // Maximum time
   const int t_pw = 100;            // Write to file every t_pw timesteps
   // Simulation Scale and Sizing
   const int N = 5;                    // Number of atoms
   const double Lx = 10.0;             // Box length in x
   const double Ly = 10.0;             // Box length in y
   const double Lz = 10.0;             // Box length in z
   // Atomic, Thermodynamic, and Molecular Conditions
   const double KbT = 1.0;             // Boltzmann constant
   const double gamma = 1.0;           // Friction coefficient
   const double length_scale = 1.0;    // Length scale for Lennard-Jones potential
   const double energy_scale = 1.0;    // Energy scale for Lennard-Jones potential
   const double rc1 = 2.5;             // Cutoff distance for Lennard-Jones potential
   const double rc2 = 2.5;             // Overlap reject thershold for random position generation
   const double T = 1.0;               // Temperature
   const double m = 1.0;               // Mass of atoms
   // Paths to write outputted data into
   std::string output_file_path_x = "./data_pos.txt";  // Data file output for positions
   std::string output_file_path_v = "./data_vel.txt";  // Data file output for velocities
/* ------------------------------------------------------------------------------------------------
   ---- Declaring Necessary Data from Input Conditions --------------------------------------------
   --------------------------------------------------------------------------------------------- */
   const int n_timesteps = (t_max - t) / dt;
   std::vector<double> timesteps(n_timesteps, 0.0);
   std::vector<double> iR(N, 0.0);
   const std::array<double, 3> L = {Lx, Ly, Lz};
   std::vector<std::vector<double>> R(3, std::vector<double>(N, 0.0));
   std::vector<std::vector<double>> V(3, std::vector<double>(N, 0.0));
/* ------------------------------------------------------------------------------------------------
   ---- Calling Initialisation Functions ----------------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   std::cout << "|| --------- Simulation and Particle Initialisation ----------- ||" << std::endl;
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   initialise_positions(R, N, L, rc2);
   initialise_velocities(V, N, L, T, m);
   pbc(R, iR, L, N);
/* ------------------------------------------------------------------------------------------------
   ---- Performing Verlet Integration -------------------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   //std::vector<std::vector<double>> force_array = force_routine(R, V, m, \
                                    gamma, KbT, dt, energy_scale, length_scale, rc1, L, N);
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   std::cout << "|| ------------- Beginning Verlet Integration ----------------- ||" << std::endl;
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   Verlet_Integration(R, V, m, N, dt, t, t_max, n_timesteps, timesteps, \
                      output_file_path_x, output_file_path_v, t_pw,\
                      gamma, KbT, energy_scale, length_scale, rc1, L, iR);

   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
   std::cout << "|| ----------- Successful Exit to End of Program -------------- ||" << std::endl;
   std::cout << "|| ------------------------------------------------------------ ||" << std::endl;
/* ------------------------------------------------------------------------------------------------
   ---- Total Wall Time ---------------------------------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> duration = end - start;
   std::cout << "Total wall time: " << duration.count() << " seconds" << std::endl;

   return 0;
}