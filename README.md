# SUPAASP

This is the repository where I have referenced the code I developed for one of my PhD grad classes (Advanced Statistical Physics). Only the MD code itself was really necessary for this class, however I realised this simple MD code would serve as an ideal 'toy-model' testing ground for exploring KOKKOS for my own learnings sake. I'm focusing on implementing a KOKKOS port for the CG-DNA (oxDNA) package within LAMMPS, though this introduces its own additional layers of abstractions and complexity - the MD code for this grad class was quite useful in providing me with an opportunity to try out a basic CPU (Host) KOKKOS implementation on a standalone package outwith LAMMPS.

For building the Vanilla C++ Serial code, this is done by running 'make' within: ./SUPAASP/src

For building the KOKKOS Port, this is done by running 'make' within: ./SUPAASP/src/cmake_build

More precisely, in my case:

> cd *ARBITRARY-KOKKOS-DOWNLOAD-DIRECTORY*/kokkos-4.3.01/ \
> mkdir BUILD && mkdir INSTALL
>> cmake -B ./BUILD/ -DKokkos_4301=ON -DCMAKE_BUILD_TYPE=Release \
>> -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=./INSTALL/ \
>> -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_ZEN3=ON -S .
>> 
> cmake --build ./BUILD/ \
> cmake --install ./BUILD/

Then, after the KOKKOS library is installed, it is linked to the project via:

> cd *ARBITRARY-PROJECT-LOCAL-DIRECTORY*/SUPAASP/src/cmake_build
>> cmake -S ./../ -B . \
>> -DKokkos_ROOT=/media/lewis/PhD/kokkos-4.3.01/INSTALL/lib/cmake/Kokkos \
>> -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release
>>
> make

