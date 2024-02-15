/* ------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   --------------------------------------------------------------------------------------------- */

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

/* --------------------------------------------------------------------------------------------- */

void prep_data_files(std::string, const int&, const int&, const std::array<double, 3>&);

void write_to_file(std::string file_write_path, std::vector<std::vector<double>>&, double, int);