/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines.h"

/* ------------------------------------------------------------------------------------------------------------------------------- */
void initialise_positions(std::vector<std::vector<double>>& R, const int& N, const std::array<double, 3>& L, const double& rc2) 
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
    printf("Initial position of particle %d: (%f, %f, %f)\n", a, xp[0], xp[1], xp[2]);
    // Check for overlap, accept pass and reject fail
    bool flag = true;
    for (int b = 0; b < N; b++) {
        double dr2 = std::pow(xp[0], 2) + std::pow(xp[1], 2) + std::pow(xp[2], 2);
        if (dr2 < rc2) {
            flag = false;
            break;
        }
    }
    if (flag) {
        R[0][a] = xp[0];
        R[1][a] = xp[1];
        R[2][a] = xp[2];
        printf("Particle %d has been initialised\n", a);
        a++;
    }
  }
  printf("Total number of particles initialised: %d\n", a);
}

/* ------------------------------------------------------------------------------------------------------------------------------- */
void initialise_velocities(std::vector<std::vector<double>>& V, const int& N, const std::array<double, 3>& L, \
                                                         const double& T, const double& m) 
{ 
  // This code section below generates random initial velocities which have the correct mean speed
  // as predicted via kinetic theory. The velocities are then scaled to the correct temperature.
  std::default_random_engine random_engine;
  for (int n = 0; n < N; n++) {
    double vnorm2 = 0.0;
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
  printf("Initial velocities have also been initialised for all %d particles\n", N);
}

/* ------------------------------------------------------------------------------------------------------------------------------- */
void pbc(std::vector<std::vector<double>>& R, std::vector<double>&iR, const std::array<double, 3>& L, const int& N) 
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

/* ------------------------------------------------------------------------------------------------------------------------------- */
std::vector<std::vector<double>> force_routine(const std::vector<std::vector<double>>& R, const std::vector<std::vector<double>>& V,\
                                               const double& m, const double& gamma, const double& KbT, const double& dt, \
                                               const double& energy_scale, const double& length_scale, const double& rc1,\
                                                std::array<double, 3>& L, const int& N) 
{
  std::vector<std::vector<double>> force_array(3, std::vector<double>(N, 0.0));

  std::vector<double> r(3, 0.0);
  std::vector<double> F_lj(3, 0.0);
  std::vector<double> eta(3, 0.0);
  std::vector<double> F_r(3, 0.0);
  std::vector<double> F_f(3, 0.0);

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
        if (std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) < rc1) {
          if (j != k) {
            F_lj[p] = -(r[p] * (4 * energy_scale * ( \
                     (std::pow(length_scale, 12) * (-12) / std::pow(std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]), 14)) \
                     + (std::pow(length_scale, 6) * (6) / std::pow(std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]), 8)))));
                     //F[j][k] = -F[k][j];
          } else {
            F_lj[p] = 0;
          }
        }
      std::random_device rd;
      std::mt19937 gen(rd());
      std::normal_distribution<double> dist(0.0, 1.0);
      eta[p] = dist(gen);
      F_r[p] = std::sqrt(2 * m * gamma * KbT) * eta[p] / std::sqrt(dt / 2);
      //F_f[p] = gamma * V[p][k];

      force_array[p][k] = ( F_lj[p]/m - gamma * V[p][k] + F_r[p]/std::sqrt(m) );
      force_array[p][j] = (-F_lj[p]/m - gamma * V[p][j] + F_r[p]/std::sqrt(m) );
      }
    }
  }
  return force_array;
}

/* ------------------------------------------------------------------------------------------------------------------------------- */
