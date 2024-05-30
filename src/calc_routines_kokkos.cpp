/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines_kokkos.h"
#include "files.h"

/* --------------------------------------------------------------------------------------------- */
/*std::vector<std::vector<double>> calc_routines_kokkos::force_routine(
  const std::vector<std::vector<double>>& R, const std::vector<std::vector<double>>& V,\
  const double& m, const double& gamma, const double& KbT, const double& dt,\
  const double& energy_scale, const double& length_scale, const double& rc1,\
  const std::array<double, 3>& L, const int& N,\
  Kokkos::Random_XorShift64_Pool<>& rand_pool)
  //Kokkos::Random_XorShift64_Pool<>::generator_type& RNG) 
{
  Kokkos::View<double*[3]> force_array_KK("force_array", N);
  //printf("Calculating forces...\n");
  Kokkos::View<double[3]> r("r",0.0);
  Kokkos::View<double[3]> F_lj("F_lj",0.0);
  Kokkos::View<double[3]> F_r("F_r",0.0);
  //Kokkos::View<double*[3]> eta(3, 0.0);

  Kokkos::parallel_for(N, KOKKOS_LAMBDA(int k) {
    for (int j = k; j < N; ++j) {
      // Minimum Image Convention
      for (int p = 0; p < 3; ++p) {
        r(p) = R[p][k] - R[p][j];
        if (r(p) > 0.5 * L[p]) {
          r(p) -= L[p];
        }
        if (r[p] < -0.5 * L[p]) {
          r(p) += L[p];
        }
      }
      double r_norm = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
      for (int p = 0; p < 3; ++p) {
        if ((r_norm < rc1) && (j != k)) {
          // Force Calculation : Lennard-Jones Potential
          F_lj(p) = -(r(p) * (4 * energy_scale * ( \
                     ((std::pow(length_scale, 12) * (-12)) / std::pow(r_norm, 14)) + \
                     ((std::pow(length_scale, 6) * (6))   /  std::pow(r_norm, 8)) )));
          //printf("F_lj(p): %f\n", F_lj(p));
          // Force Calculation: Langevin Dynamics
          //printf("Calculating random number...\n");
          auto RNG = rand_pool.get_state();
          double eta = RNG.drand(0.0, 1.0);
          rand_pool.free_state(RNG);
          eta = 0.0;
          F_r(p) = std::sqrt(2 * m * gamma * KbT) * eta / std::sqrt(dt / 2);
          //printf("F_r(p): %f\n", F_r(p));
          //printf("F_r: Done\n");
          // Force Calculation: Total Force
          force_array_KK(k,p) = ( F_lj(p)/m - gamma * V[p][k] + F_r(p)/std::sqrt(m) );
          force_array_KK(j,p) = (-F_lj(p)/m - gamma * V[p][j] + F_r(p)/std::sqrt(m) );
          //printf("force_array_KK(p,k): %f\n", force_array_KK(k,p));
          //printf("force_array_KK(p,j): %f\n", force_array_KK(k,j));
          //printf("force_array_KK increment: Done\n");
        }
      }
      //printf("End of inner loop\n");
    }
    //printf("End of outer loop\n");
  });
  //rand_pool.free_state(RNG);
  //return force_array_KK;
  //printf("End of parallel_for\n");

  return force_array_KK;
}*/

/* --------------------------------------------------------------------------------------------- */
void calc_routines_kokkos::Verlet_Integration(
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
  //std::vector<std::vector<double>> temp_force_routine(3, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> temp_force_routine(3, std::vector<double>(N, 0.0));
  std::default_random_engine random_engine(seed);
  std::normal_distribution<double> dist(0.0, 1.0);

  //Kokkos::View<double*[3]> temp_force_routine("temp_force_routine", N);
  //printf("Initialising random number generator pool...\n");
  //Kokkos::Random_XorShift64_Pool<> rand_pool(seed);
  //Kokkos::rand<Kokkos::Random_XorShift64_Pool<>::generator_type, double> rand_gen(rand_pool);
  //Kokkos::Random_XorShift64_Pool<>::generator_type RNG = rand_pool.get_state();
  //printf("Beginning Verlet Integration...\n");
  for (int i = 0; i <= t_n; i++) {
    temp_force_routine = force_routine(R, V, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, \
                                       //rand_pool);
                                       random_engine, dist);

    Kokkos::parallel_for(N, KOKKOS_LAMBDA(int n) {
      for (int p = 0; p < 3; p++) {
        V_i_half[p][n] = V[p][n] + temp_force_routine[p][n] * dt / 2;
        R[p][n] += V_i_half[p][n] * dt;
      }
    });
    
    temp_force_routine = force_routine(R, V_i_half, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, \
                                       //rand_pool);
                                       random_engine, dist);

    Kokkos::parallel_for(N, KOKKOS_LAMBDA(int n) {
      for (int p = 0; p < 3; p++) {
        V[p][n] = V_i_half[p][n] + temp_force_routine[p][n] * dt / 2;
      }
    });

    pbc(R, iR, L, N);
    if ( (i % t_pw) == 0) {
      printf("\rWriting to data file at iteration: %d, time: %f", i, t_i);
      write_to_file(filename_x, R, t_i, i); // store (appending) the position data to file 
      write_to_file(filename_v, V, t_i, i); // store (appending) the velocity data to file
    }
    timesteps[i] = t_i; // stores time values
    t_i += dt;
  }
}