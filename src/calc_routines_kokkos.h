/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

class calc_routines_kokkos : public calc_routines {
public:

   calc_routines_kokkos(){};  // Empty constructor
   ~calc_routines_kokkos(){}; // Empty destructor

   void write_to_file_KK(std::string file_write_path, const Kokkos::View<double**>,\
                     double, int, const int&);

   Kokkos::View<double**> force_routine(const Kokkos::View<double**>&,\
                     const Kokkos::View<double**>&,\
                     const double&, const double&, const double&,\
                     const double&, const double&, const double&,\
                     const double&, const std::array<double, 3>&, const int&,\
                     Kokkos::Random_XorShift64_Pool<>&);

   void Verlet_Integration(std::vector<std::vector<double>>, std::vector<std::vector<double>>,\
                     const double&, const int&, const double& dt, const double& t,\
                     const double&, const int&, std::vector<double>&,\
                     const std::string&, const std::string&, const int&,\
                     const double&, const double&, const double&, const double&,\
                     const double&, const std::array<double, 3>&, std::vector<double>&,\
                     const int&) override;
};
