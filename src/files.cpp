/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "files.h"
#include <iomanip> // Include the <iomanip> library for std::setw()
                   // setw() ensures column alignment - looks much tidier.

/* --------------------------------------------------------------------------------------------- */
void prep_data_files(std::string file_write_path, const int& N, const int& n_timesteps, \
                     const std::array<double, 3>& L) {
   std::ofstream file_out;
   file_out.open(file_write_path, std::ios::trunc); // Opens the file in truncation mode
   if (file_out.fail()) {
      std::cout << "ERROR: Could not open or write to file!\n" << std::endl;
   }
   file_out << "Number of Particles: " << N << std::endl;
   file_out << "Number of Timesteps: " << n_timesteps << std::endl;
   file_out << "Box Lengths (x, y, z): " << std::setw(2) \
            << L[0] << " , " << L[1] << " , " << L[2] << std::endl;
   file_out << "\nTimestep results below:\n" <<std::endl;
   file_out.close();
   std::cout << "File \"" <<  \
          file_write_path << "\" has been created or wiped, now prepared!" << std::endl;
}

/* --------------------------------------------------------------------------------------------- */
void write_to_file(std::string file_write_path, std::vector<std::vector<double>>& \
                  data_vec_3D, double t, int it) {

   std::ofstream file_out;
   file_out.open(file_write_path, std::ios::app); // Opens the file in append mode
   if (file_out.fail()) {
      std::cout << "ERROR: Could not open or write to file!\n" << std::endl;
   } else {
      //file_out << "\nData Points (x \\t \\t y \\t \\t z):\n" <<std::endl;
      int n = data_vec_3D[1].size();
      file_out << "Time: " << t << std::endl;
      file_out << "Iteration Step: " << it << std::endl;
      for (int i=0; i<n; i++) {
         file_out << std::setw(10) << data_vec_3D[0][i] << \
             "\t" << std::setw(10) << data_vec_3D[1][i] << \
             "\t" << std::setw(10) << data_vec_3D[2][i] << std::endl;
      }
   file_out.close();
   std::cout << "Successful data append to file \"" << file_write_path << "\"!" << std::endl;
   }
}