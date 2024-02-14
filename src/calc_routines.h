/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <random>

void pbc(std::vector<std::vector<double>>&, std::vector<double>&, \
                                               const std::array<double, 3>&, const int&);

std::vector<std::vector<double>> force_routine(const std::vector<std::vector<double>>&,\
                                               const std::vector<std::vector<double>>&,\
                                               const double&, const double&, const double&,\
                                               const double&, const double&, const double&,\
                                               const double&, std::array<double, 3>&,\
                                               const int&);

void initialise_positions(std::vector<std::vector<double>>&,\
                          const int&, const std::array<double, 3>&, const double&);

void initialise_velocities(std::vector<std::vector<double>>&,\
                           const int&, const std::array<double, 3>&, \
                           const double&, const double&);