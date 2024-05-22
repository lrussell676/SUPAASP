/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines_kokkos.h"
#include "files.h"

/* --------------------------------------------------------------------------------------------- */
template<class Device>
std::vector<std::vector<double>> calc_routines_kokkos<Device>::force_routine(
  const std::vector<std::vector<double>>& R, const std::vector<std::vector<double>>& V,\
  const double& m, const double& gamma, const double& KbT, const double& dt,\
  const double& energy_scale, const double& length_scale, const double& rc1,\
  const std::array<double, 3>& L, const int& N,\
  std::default_random_engine& random_engine, std::normal_distribution<double>& dist) 
{
  std::vector<std::vector<double>> force_array(3, std::vector<double>(N, 0.0));

  std::vector<double> r(3, 0.0);
  std::vector<double> F_lj(3, 0.0);
  std::vector<double> F_r(3, 0.0);
  std::vector<double> eta(3, 0.0);
  double r_norm;

  printf("USING KOKKOS\n");

  for (int k = 0; k < N; ++k) {
    for (int j = k; j < N; ++j) {
      for (int p = 0; p < 3; ++p) {
        // Minimum Image Convention
        r[p] = R[p][k] - R[p][j];
        if (r[p] > 0.5 * L[p]) {
          r[p] -= L[p];
        }
        if (r[p] < -0.5 * L[p]) {
          r[p] += L[p];
        }
        // Force Calculation : Lennard-Jones Potential
        r_norm = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        if ( (r_norm < rc1) && (j != k) ) {
          F_lj[p] = -(r[p] * (4 * energy_scale * ( \
                     ((std::pow(length_scale, 12) * (-12)) / std::pow(r_norm, 14)) + \
                     ((std::pow(length_scale, 6) * (6))   /  std::pow(r_norm, 8)) )));
        } else F_lj[p] = 0.0;
        // Force Calculation: Langevin Dynamics
        eta[p] = dist(random_engine);
        F_r[p] = std::sqrt(2 * m * gamma * KbT) * eta[p] / std::sqrt(dt / 2);
        // Force Calculation: Total Force
        force_array[p][k] = ( F_lj[p]/m - gamma * V[p][k] + F_r[p]/std::sqrt(m) );
        force_array[p][j] = (-F_lj[p]/m - gamma * V[p][j] + F_r[p]/std::sqrt(m) );
      }
    }
  }
  return force_array;
}

/* --------------------------------------------------------------------------------------------- *
void calc_routines::Verlet_Integration(
  std::vector<std::vector<double>> R, std::vector<std::vector<double>> V, const double& m,\
  const int& N, const double& dt, const double& t, const double& t_max, const int& t_n,\
  std::vector<double>& timesteps, const std::string& filename_x, const std::string& filename_v,\
  const int& t_pw, const double& gamma, const double& KbT, const double& energy_scale,\
  const double& length_scale, const double& rc1,\
  const std::array<double, 3>& L, std::vector<double>& iR, const int& seed) 
{
  double t_i = t;
  std::vector<double> time;
  std::vector<std::vector<double>> V_i_half(3, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> temp_force_routine(3, std::vector<double>(N, 0.0));
  std::default_random_engine random_engine(seed);
  std::normal_distribution<double> dist(0.0, 1.0);
  
  for (int i = 0; i <= t_n; i++) {
    temp_force_routine = force_routine(R, V, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, \
                                       random_engine, dist);
    for (int n = 0; n < N; n++) {
      for (int p = 0; p < 3; p++) {
        V_i_half[p][n] = V[p][n] + temp_force_routine[p][n] * dt / 2;
        R[p][n] += V_i_half[p][n] * dt;
      }
    }
    temp_force_routine = force_routine(R, V_i_half, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, \
                                       random_engine, dist);
    for (int n = 0; n < N; n++) {
      for (int p = 0; p < 3; p++) {
        V[p][n] = V_i_half[p][n] + temp_force_routine[p][n] * dt / 2;
      }
    }
    pbc(R, iR, L, N);
    if ( (i % t_pw) == 0) {
      printf("\rWriting to data file at iteration: %d, time: %f", i, t_i);
      write_to_file(filename_x, R, t_i, i); // store (appending) the position data to file 
      write_to_file(filename_v, V, t_i, i); // store (appending) the velocity data to file
    }
    timesteps[i] = t_i; // stores time values
    t_i += dt;
  }
}*/