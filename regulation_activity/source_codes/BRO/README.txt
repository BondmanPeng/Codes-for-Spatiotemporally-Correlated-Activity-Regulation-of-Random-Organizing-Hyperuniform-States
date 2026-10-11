Build from this directory:
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build -j2
  ctest --test-dir build --output-on-failure
Requires C++17 and OpenMP (BRO); GCC 12.3 was used to check this package.
The simulation default h0 is 2; a supplied positional h0 still overrides it.
Exact circular density: minimum-image particle distances <= R_A, divided by pi*R_A^2.
Cells accelerate candidate lookup only; they do not approximate disk membership.

Figure data are packaged as CSVs in ../../available_data/BRO.
Plotting scripts are not currently included in this package.
See ../../README.txt for dataset details and original source locations.

2026-10-09: rho_A and rho_B are now number densities: N=floor(rho*L^2).
At L=300, rho_A=0.8 and rho_B=0.6 give N_A=72000, N_B=54000.
Exact circular sensing and the simulation default h0=2 are retained.
The packaged BRO figure tables now use rho_A=0.8, rho_B=0.7 and explicit h0=1.5
with number-density conventions. See ../../README.txt for selections and provenance.
