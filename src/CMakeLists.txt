# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.21.1)
set(CMAKE_CXX_STANDARD_REQUIRED 20)

# Set the project name
project(SUPAASP_KOKKOS)

# Add the source files
file(GLOB SOURCES "*.cpp" "*.h")

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../src)

# Link to the installed package
find_package(Kokkos REQUIRED)
if(Kokkos_FOUND)
    message("-- Kokkos found")
    target_link_libraries(${PROJECT_NAME}  Kokkos::kokkos)
    add_compile_definitions(KOKKOS_ENABLED)
else()
  message(FATAL_ERROR "Kokkos not found")
endif()


#macro(check_package package)
#  find_package(${package})
#  if(${package}_FOUND)
#    message("Package ${package} found")
#  else()
#    message("Package ${package} not found")
#  endif()
#endmacro()
