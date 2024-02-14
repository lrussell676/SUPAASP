/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include "files.h"

void write_to_file(std::string file_write_path, std::vector<std::vector<double>>& \
                  data_vec_3D) {

   std::ofstream file_out;

   file_out.open(file_write_path);

   if (file_out.fail()) {
      std::cout << "ERROR: Cound not open or write to file.\n" << std::endl;
   } else {
      file_out << "\nData Points (x \\t \\t y \\t \\t z):\n" <<std::endl;
      int n = data_vec_3D[1].size();
      for (int i=0; i<n; i++) {
         file_out << \
         data_vec_3D[0][i] << "\t\t" << data_vec_3D[1][i] << "\t\t" << data_vec_3D[2][i] \
         << std::endl;
      }
   file_out.close();
   std::cout << "Successful write to file!\n" << std::endl;
   }
   
}