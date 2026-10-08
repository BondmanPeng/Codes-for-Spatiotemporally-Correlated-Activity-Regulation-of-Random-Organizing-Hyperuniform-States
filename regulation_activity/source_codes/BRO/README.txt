Build from this directory:
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build -j2
  ctest --test-dir build --output-on-failure
Requires C++17 and OpenMP (BRO); GCC 12.3 was used to check this package.
The simulation default h0 is 2; a supplied positional h0 still overrides it.
Exact circular density: minimum-image particle distances <= R_A, divided by pi*R_A^2.
Cells accelerate candidate lookup only; they do not approximate disk membership.

Plotting (requires Python 3, NumPy and Matplotlib):
  python3 plotting/plot_combined_sbb_c_h.py
  python3 plotting/plot_saa_c_sweep.py
  python3 plotting/plot_delta_curves.py

Plotters use paths relative to this package and write PNG/PDF alongside their
input tables in available_data. No simulation or scheduler submission is run.
See ../../README.txt for dataset and density-convention details.
