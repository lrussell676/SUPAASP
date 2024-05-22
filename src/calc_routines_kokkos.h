/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "calc_routines.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Timer.hpp>
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_UnorderedMap.hpp>

template<class Device>
class calc_routines_kokkos : public calc_routines {
public:

   calc_routines_kokkos(){};  // Empty constructor
   ~calc_routines_kokkos(){}; // Empty destructor

   //void pbc(std::vector<std::vector<double>>&, std::vector<double>&,\
   //       const std::array<double, 3>&, const int&);

   KOKKOS_INLINE_FUNCTION
   std::vector<std::vector<double>> force_routine(const std::vector<std::vector<double>>&,\
                                       const std::vector<std::vector<double>>&,\
                                       const double&, const double&, const double&,\
                                       const double&, const double&, const double&,\
                                       const double&, const std::array<double, 3>&,\
                                       const int&, std::default_random_engine&,\
                                       std::normal_distribution<double>&) override;

   /*void initialise_positions(std::vector<std::vector<double>>&,\
                       const int&, const std::array<double, 3>&, const double&);

   void initialise_velocities(std::vector<std::vector<double>>&,\
                        const int&, const std::array<double, 3>&,\
                        const double&, const double&, const int&);

   void Verlet_Integration(std::vector<std::vector<double>>, std::vector<std::vector<double>>,\
                     const double&, const int&, const double& dt, const double& t,\
                     const double&, const int&, std::vector<double>&,\
                     const std::string&, const std::string&, const int&,\
                     const double&, const double&, const double&, const double&,\
                     const double&, const std::array<double, 3>&, std::vector<double>&,\
                     const int&);*/
};
