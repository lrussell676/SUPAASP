/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines_kokkos.h"
#include "files.h"

/* --------------------------------------------------------------------------------------------- */
void calc_routines_kokkos::write_to_file_KK(std::string file_write_path, Kokkos::View<double**> \
                  data_vec_3D, double t, int it, const int& n) {

   std::ofstream file_out;
   file_out.open(file_write_path, std::ios::app); // Opens the file in append mode
   if (file_out.fail()) {
      std::cout << "ERROR: Could not open or write to file!\n" << std::endl;
   } else {
      file_out << "Time: " << t << std::endl;
      file_out << "Iteration Step: " << it << std::endl;
      for (int i=0; i<n; i++) {
         file_out << std::setw(10) << data_vec_3D(0,i) << \
             "\t" << std::setw(10) << data_vec_3D(1,i) << \
             "\t" << std::setw(10) << data_vec_3D(2,i) << std::endl;
      }
   file_out.close();
   }
}

/* --------------------------------------------------------------------------------------------- */
Kokkos::View<double**> calc_routines_kokkos::force_routine(
  const Kokkos::View<double**>& R, const Kokkos::View<double**>& V,\
  const double& m, const double& gamma, const double& KbT, const double& dt,\
  const double& energy_scale, const double& length_scale, const double& rc1,\
  const std::array<double, 3>& L, const int& N,\
  Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
  Kokkos::View<double**> force_array_KK("force_array", R.extent(0), R.extent(1));
  typedef typename Kokkos::Random_XorShift64_Pool<>::generator_type rand_type;
  // Parallel Force Routine
  Kokkos::parallel_for(N, KOKKOS_LAMBDA(int k) {
    double F_lj, F_r;
    Kokkos::View<double[3]> r("r",0.0);
    Kokkos::View<double[3]> eta("eta",0.0);
    for (int j = k; j < N; ++j) {
      // Minimum Image Convention
      for (int p = 0; p < 3; ++p) {
        r(p) = R(p,k) - R(p,j);
        if (r(p) > 0.5 * L[p]) {
          r(p) -= L[p];
        }
        if (r(p) < -0.5 * L[p]) {
          r(p) += L[p];
        }
      }
      double r_norm = std::sqrt(r(0) * r(0) + r(1) * r(1) + r(2) * r(2));
      rand_type RNG = rand_pool.get_state();
      eta(0) = RNG.normal(0.0,1.0);
      eta(1) = RNG.normal(0.0,1.0);
      eta(2) = RNG.normal(0.0,1.0);
      rand_pool.free_state(RNG);
      for (int p = 0; p < 3; ++p) {
        // Force Calculation : Lennard-Jones Potential
        if ((r_norm < rc1) && (j != k)) {
          F_lj = -(r(p) * (4 * energy_scale * ( \
                     ((std::pow(length_scale, 12) * (-12)) / std::pow(r_norm, 14)) + \
                     ((std::pow(length_scale, 6) * (6))   /  std::pow(r_norm, 8)) )));
        } else F_lj = 0.0;
        // Force Calculation: Langevin Dynamics
        F_r = std::sqrt(2 * m * gamma * KbT) * eta(p) / std::sqrt(dt / 2);
        // Force Calculation: Total Force
        force_array_KK(p,k) = ( F_lj/m - gamma * V(p,k) + F_r/std::sqrt(m) );
        force_array_KK(p,j) = (-F_lj/m - gamma * V(p,j) + F_r/std::sqrt(m) );
      }
    }
  });
  // End of Parallel Force Routine
  return force_array_KK;
}

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
  // Initialising Kokkos Views
  Kokkos::View<double**> V_i_half_KK("V_i_half_KK", V.size(), V[0].size());
  Kokkos::View<double**> temp_force_routine("temp_force_routine", V.size(), V[0].size());
  Kokkos::View<double**> R_KK("R_KK", R.size(), R[0].size());
  Kokkos::View<double**> V_KK("V_KK", V.size(), V[0].size());
  // Copying Vanilla C++ std::vectors to Kokkos Views
  for (int p = 0; p < 3; p++) {
    for (int n = 0; n < N; n++) {
      R_KK(p,n) = R[p][n];
      V_KK(p,n) = V[p][n];
    }
  }
  // Defining Kokkos RNG Pool
  Kokkos::Random_XorShift64_Pool<> rand_pool(seed);
  // Verlet Integration Loop
  for (int i = 0; i <= t_n; i++) {
    // 1st Force Routine Call
    temp_force_routine = force_routine(R_KK, V_KK, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, \
                                       rand_pool);
    // Parallel Verlet 1st Step
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(int n) {
      for (int p = 0; p < 3; p++) {
        V_i_half_KK(p,n) = V_KK(p,n) + temp_force_routine(p,n) * dt / 2;
        R_KK(p,n) += V_i_half_KK(p,n) * dt;
      }
    });
    // 2nd Force Routine Call
    temp_force_routine = force_routine(R_KK, V_i_half_KK, m, gamma, KbT, dt, \
                                       energy_scale, length_scale, rc1, L, N, \
                                       rand_pool);
    // Parallel Verlet 2nd Step
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(int n) {
      for (int p = 0; p < 3; p++) {
        V_KK(p,n) = V_i_half_KK(p,n) + temp_force_routine(p,n) * dt / 2;
      }
    });
    // Periodic Boundary Conditions
    for (int n = 0; n < N; ++n) {
      for (int p = 0; p < 3; ++p) {
        if (R_KK(p,n) < -0.5 * L[p]) {
          R_KK(p,n) += L[p];
          iR[n] -= 1;
        }
        if (R_KK(p,n) > 0.5 * L[p]) {
          R_KK(p,n) -= L[p];
          iR[n] += 1;
        }
      }
    }
    // Writing to Data Files
    if ( (i % t_pw) == 0) {
      printf("\rWriting to data file at iteration: %d, time: %f", i, t_i);
      write_to_file_KK(filename_x, R_KK, t_i, i, N); // store (appending) the position data to file 
      write_to_file_KK(filename_v, V_KK, t_i, i, N); // store (appending) the velocity data to file
    }
    timesteps[i] = t_i; // stores time values
    t_i += dt;
  }
}