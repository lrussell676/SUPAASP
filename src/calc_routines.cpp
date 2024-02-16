/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines.h"
#include "files.h"

/* --------------------------------------------------------------------------------------------- */
void initialise_positions(std::vector<std::vector<double>>& R, const int& N, \
                          const std::array<double, 3>& L, const double& rc2) 
{
  std::default_random_engine random_engine;
  int a = 0;
  while (a < N) {
    // Initialise random positions with Cartesian components in (-1, +1)
    std::array<double, 3> xp = {
        std::uniform_real_distribution<double>(-1.0, 1.0)(random_engine),
        std::uniform_real_distribution<double>(-1.0, 1.0)(random_engine),
        std::uniform_real_distribution<double>(-1.0, 1.0)(random_engine)
    };
    // Scale random coordinates to box size
    xp[0] *= 0.5 * L[0];
    xp[1] *= 0.5 * L[1];
    xp[2] *= 0.5 * L[2];
    // Check for overlap, accept pass and reject fail
    bool flag = true;
    for (int b = 0; b < N; b++) {
        double dr2 = std::sqrt(std::pow(xp[0]-R[0][b],2) + \
                               std::pow(xp[1]-R[1][b], 2) + std::pow(xp[2]-R[2][b], 2));
        if (dr2 < rc2) {
            flag = false;
            printf("Overlap detected, rejecting and reinitialising particle %d\n", a);
            break;
        }
    }
    if (flag) {
        R[0][a] = xp[0];
        R[1][a] = xp[1];
        R[2][a] = xp[2];
        a++;
    }
  }
  printf("Total number of particles initialised: %d\n", a);
}

/* --------------------------------------------------------------------------------------------- */
void initialise_velocities(std::vector<std::vector<double>>& V, const int& N, \
                           const std::array<double, 3>& L, const double& T, const double& m, \
                           const int& seed) 
{ 
  // This code section below generates random initial velocities which have the correct mean speed
  // as predicted via kinetic theory. The velocities are then scaled to the correct temperature.
  std::default_random_engine random_engine(seed);
  for (int n = 0; n < N; n++) {
    double vnorm2;
    for (int i = 0; i < 3; i++) {
        V[i][n] = std::uniform_real_distribution<double>(-1.0, 1.0)(random_engine);
        vnorm2 += std::pow(V[i][n], 2);
    }
    double vnorm = std::sqrt(vnorm2);
    for (int i = 0; i < 3; i++) {
        V[i][n] /= vnorm;
        V[i][n] *= std::sqrt(3 * T / m);
    }
  }
  printf("Initial velocities have also been initialised for all particles.\n");
}

/* --------------------------------------------------------------------------------------------- */
void pbc(std::vector<std::vector<double>>& R, std::vector<double>&iR, \
         const std::array<double, 3>& L, const int& N) 
{
  for (int n = 0; n < N; ++n) {
    for (int i = 0; i < 3; ++i) {
      if (R[i][n] < -0.5 * L[i]) {
        R[i][n] += L[i];
        iR[n] -= 1;
      }
      if (R[i][n] > 0.5 * L[i]) {
        R[i][n] -= L[i];
        iR[n] += 1;
      }
    }
  }
}

/* --------------------------------------------------------------------------------------------- */
std::vector<std::vector<double>> force_routine(
  const std::vector<std::vector<double>>& R, const std::vector<std::vector<double>>& V,\
  const double& m, const double& gamma, const double& KbT, const double& dt, \
  const double& energy_scale, const double& length_scale, const double& rc1,\
  const std::array<double, 3>& L, const int& N, const int& seed) 
{
  std::vector<std::vector<double>> force_array(3, std::vector<double>(N, 0.0));

  std::vector<double> r(3, 0.0);
  std::vector<double> F_lj(3, 0.0);
  std::vector<double> eta(3, 0.0);
  std::vector<double> F_r(3, 0.0);
  std::vector<double> F_f(3, 0.0);
  double r_norm;

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
        //printf("Distance between particles %d and %d: %f\n", k, j, r_norm);
        if ( (r_norm < rc1) && (j != k) ) {
          F_lj[p] = -(r[p] * (4 * energy_scale * ( \
                     ((std::pow(length_scale, 12) * (-12)) / std::pow(r_norm, 14)) + \
                     ((std::pow(length_scale, 6) * (6))   /  std::pow(r_norm, 8)) )));
        } else F_lj[p] = 0.0;
        //printf("Lennard-Jones force component %d for particle %d: %f\n", p, k, F_lj[p]);
        std::default_random_engine random_engine(seed);
        std::normal_distribution<double> dist(0.0, 1.0);
        eta[p] = dist(random_engine);
        F_r[p] = std::sqrt(2 * m * gamma * KbT) * eta[p] / std::sqrt(dt / 2);
        //F_f[p] = gamma * V[p][k];

        force_array[p][k] = ( F_lj[p]/m - gamma * V[p][k] + F_r[p]/std::sqrt(m) );
        force_array[p][j] = (-F_lj[p]/m - gamma * V[p][j] + F_r[p]/std::sqrt(m) );
        //printf("Force component %d for particle %d: %f\n", p, k, force_array[p][k]);
        //printf("Force component %d for particle %d: %f\n", p, j, force_array[p][j]);
      }
    }
  }
  return force_array;
}

/* --------------------------------------------------------------------------------------------- */
void Verlet_Integration(
  std::vector<std::vector<double>>& R, std::vector<std::vector<double>>& V, const double& m, \
  const int& N, const double& dt, const double& t, const double& t_max, const int& t_n, \
  std::vector<double>& timesteps, const std::string& filename_x, const std::string& filename_v, \
  const int& t_pw, const double& gamma, const double& KbT, const double& energy_scale, \
  const double& length_scale, const double& rc1,\ 
  const std::array<double, 3>& L, std::vector<double>& iR, const int& seed) 
{
  double t_i = t;
  std::vector<double> time;
  std::vector<std::vector<double>> V_i_half(3, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> temp_force_routine(3, std::vector<double>(N, 0.0));
  int i = 0;
  for (int i = 0; i <= t_n; i++) {
    temp_force_routine = force_routine(R, V, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, seed);
    for (int n = 0; n < N; n++) {
      for (int p = 0; p < 3; p++) {
        V_i_half[p][n] = V[p][n] + temp_force_routine[p][n] * dt / 2;
        R[p][n] += V_i_half[p][n] * dt;
      }
    }
    temp_force_routine = force_routine(R, V_i_half, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, seed);
    for (int n = 0; n < N; n++) {
      for (int p = 0; p < 3; p++) {
        V[p][n] = V_i_half[p][n] + temp_force_routine[p][n] * dt / 2;
      }
    }
    pbc(R, iR, L, N);
    if ( (i % t_pw) == 0) {
      printf("Writing to data file at iteration: %d, time: %f\n", t_i, i);
      write_to_file(filename_x, R, t_i, i); // store (appending) the position data to file 
      write_to_file(filename_v, V, t_i, i); // store (appending) the velocity data to file
    }
    timesteps[i] = t_i; // stores time values
    t_i += dt;
  }
}