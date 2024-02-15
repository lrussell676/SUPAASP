/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "files.h"
#include <iomanip> // Include the <iomanip> library for std::setw()
                   // setw() ensures column alignment - looks much tidier.

/* --------------------------------------------------------------------------------------------- */
void clear_data_files(std::string file_write_path) {
   std::ofstream file_out;
   file_out.open(file_write_path, std::ios::trunc); // Opens the file in truncation mode
   file_out.close();
   std::cout << "File \"" << file_write_path << "\" has been created/cleared!" << std::endl;
}

/* --------------------------------------------------------------------------------------------- */
void write_to_file(std::string file_write_path, std::vector<std::vector<double>>& \
                  data_vec_3D) {

   std::ofstream file_out;

   file_out.open(file_write_path, std::ios::app); // Opens the file in append mode

   if (file_out.fail()) {
      std::cout << "ERROR: Could not open or write to file.\n" << std::endl;
   } else {
      //file_out << "\nData Points (x \\t \\t y \\t \\t z):\n" <<std::endl;
      int n = data_vec_3D[1].size();
      for (int i=0; i<n; i++) {
         file_out << std::setw(10) << data_vec_3D[0][i] << \
             "\t" << std::setw(10) << data_vec_3D[1][i] << \
             "\t" << std::setw(10) << data_vec_3D[2][i] << std::endl;
      }
   file_out.close();
   std::cout << "Successful data append to file \"" << file_write_path << "\"!" << std::endl;
   }
}