/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <random>

class calc_routines {
public:

   calc_routines(){};  // Empty constructor
   ~calc_routines(){}; // Empty destructor

   void initialise_positions(std::vector<std::vector<double>>&,\
                       const int&, const std::array<double, 3>&,\
                       const double&, const int&);

   void initialise_velocities(std::vector<std::vector<double>>&,\
                        const int&, const std::array<double, 3>&,\
                        const double&, const double&, const int&);

   void pbc(std::vector<std::vector<double>>&, std::vector<double>&,\
          const std::array<double, 3>&, const int&);

   virtual std::vector<std::vector<double>> force_routine(const std::vector<std::vector<double>>&,\
                                       const std::vector<std::vector<double>>&,\
                                       const double&, const double&, const double&,\
                                       const double&, const double&, const double&,\
                                       const double&, const std::array<double, 3>&,\
                                       const int&, std::default_random_engine&,\
                                       std::normal_distribution<double>&);

   virtual void Verlet_Integration(std::vector<std::vector<double>>, std::vector<std::vector<double>>,\
                     const double&, const int&, const double& dt, const double& t,\
                     const double&, const int&, std::vector<double>&,\
                     const std::string&, const std::string&, const int&,\
                     const double&, const double&, const double&, const double&,\
                     const double&, const std::array<double, 3>&, std::vector<double>&,\
                     const int&);
};
