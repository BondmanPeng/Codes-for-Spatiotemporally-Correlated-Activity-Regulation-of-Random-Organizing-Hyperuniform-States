#include "randorg_binary_randA.h"

#include <limits>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
  double L = 300.0;
  double radius_A = 0.5, radius_B = 0.5;
  double rho_A = (argc > 1) ? stod(argv[1]) : 0.60;
  double rho_B = (argc > 2) ? stod(argv[2]) : 0.60;
  double h = (argc > 3) ? stod(argv[3]) : 0.50;
  double h0 = (argc > 4) ? stod(argv[4]) : 2.0;
  double c = (argc > 5) ? stod(argv[5]) : 0.0;
  double R_A = (argc > 6) ? stod(argv[6]) : 5.00;
  int steps = (argc > 7) ? stoi(argv[7]) : 100000;
  int n_snap = (argc > 8) ? stoi(argv[8]) : 500;
  int tau_ratio = (argc > 9) ? stoi(argv[9]) : 1;
  double epsilon_A = (argc > 10) ? stod(argv[10]) : 1.0;

  const int tail_count = steps - steps / 2;
  if (tau_ratio > 0 && steps > numeric_limits<int>::max() / tau_ratio)
    throw invalid_argument("steps * tau_ratio is too large");

  ofstream dump("configuration.lammpstrj");
  if (!dump) throw runtime_error("could not open configuration.lammpstrj");

  RandOrgBinary sys;
  sys.initialize_system(L, radius_A, radius_B, rho_A, rho_B, h, h0, R_A,
                        c, tau_ratio, epsilon_A);
  printf("  mean A number density = %.4f\n", sys.mean_A_number_density);
  printf("  tau_B/tau_A = %d; total A ticks = %lld\n", tau_ratio,
         static_cast<long long>(steps) * tau_ratio);
  printf("  dumping %d snapshots over the final %d B cycles\n", n_snap,
         tail_count);

  sys.run(steps, dump, n_snap,             true);
  dump.close();
  sys.cleanup();

  printf("\nSaved configuration.lammpstrj. Structure factor: run "
         "compute_sq_species\n");
  return 0;
}
