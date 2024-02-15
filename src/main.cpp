/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines.h"
#include "files.h"

int main() {
/* ------------------------------------------------------------------------------------------------
   ---- ### DEFINE SIMULATION CONDITIONS HERE ### -------------------------------------------------
   --------------------------------------------------------------------------------------------- */
   // Timescales
   const double dt = 0.0005;           // Timestep size
   const double t = 0.0;               // Initial time
   const double t_max = 10.0;          // Maximum time
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
/* ------------------------------------------------------------------------------------------------
   ---- Declaring Necessary Vectors and Arrays ----------------------------------------------------
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
   initialise_positions(R, N, L, rc2);
   initialise_velocities(V, N, L, T, m);
   std::string output_file_path_1 = "./initial_pos.txt";
   std::string output_file_path_2 = "./initial_vel.txt";
   //write_to_file(output_file_path_1, R);
   //write_to_file(output_file_path_2, V);
   pbc(R, iR, L, N);
   std::vector<std::vector<double>> force_array = force_routine(R, V, m, \
                                    gamma, KbT, dt, energy_scale, length_scale, rc1, L, N);
   //write("filename_x.txt", "filename_v.txt", R, V);
   Verlet_Integration(R, V, m, N, dt, t, t_max, n_timesteps, timesteps, \
                      output_file_path_1, output_file_path_2,\
                      gamma, KbT, energy_scale, length_scale, rc1, L, iR);

   std::cout << "Hello, World!" << std::endl;

   return 0;
}