#include "randorg_binary_randA.h"

#include <limits>
#include <string>

using namespace std;

namespace {
bool parse_int_strict(const char *text, int &value) {
  size_t consumed = 0;
  try {
    const long long parsed = stoll(text, &consumed);
    if (consumed != string(text).size() ||
        parsed < numeric_limits<int>::min() ||
        parsed > numeric_limits<int>::max())
      return false;
    value = static_cast<int>(parsed);
    return true;
  } catch (const exception &) {
    return false;
  }
}
}

int main(int argc, char *argv[]) {
  double L = 300.0;
  double radius_A = 0.5, radius_B = 0.5;
  double rho_A = (argc > 1) ? stod(argv[1]) : 0.60;
  double rho_B = (argc > 2) ? stod(argv[2]) : 0.60;
  double h = (argc > 3) ? stod(argv[3]) : 0.50;
  double h0 = (argc > 4) ? stod(argv[4]) : 2.0;
  double c = (argc > 5) ? stod(argv[5]) : 0.0;
  double R_A = (argc > 6) ? stod(argv[6]) : 5.00;
  int steps = 10000;
  int n_snap = 100;
  int tau_ratio = 1;
  double epsilon_A = (argc > 10) ? stod(argv[10]) : 1.0;

  if (argc > 7 && !parse_int_strict(argv[7], steps)) {
    cerr << "steps must be an integer\n";
    return 2;
  }
  if (argc > 8 && !parse_int_strict(argv[8], n_snap)) {
    cerr << "n_snap must be an integer\n";
    return 2;
  }
  if (argc > 9 && !parse_int_strict(argv[9], tau_ratio)) {
    cerr << "tau_ratio must be an integer >= 0\n";
    return 2;
  }
  if (!isfinite(c) || c < -1.0 || c > 0.0) {
    cerr << "c must satisfy -1 <= c <= 0\n";
    return 2;
  }
  if (!isfinite(epsilon_A) || epsilon_A < 0.0) {
    cerr << "epsilon_A must be finite and >= 0\n";
    return 2;
  }

  const int tail_count = steps > 0 ? steps - steps / 2 : 0;
  if (steps < 1) {
    cerr << "steps must be >= 1\n";
    return 2;
  }
  if (tau_ratio < 0) {
    cerr << "tau_ratio must be >= 0\n";
    return 2;
  }
  if (n_snap < 1 || n_snap > tail_count) {
    cerr << "n_snap must satisfy 1 <= n_snap <= steps - steps/2\n";
    return 2;
  }
  if (tau_ratio > 0 && steps > numeric_limits<int>::max() / tau_ratio) {
    cerr << "steps * tau_ratio is too large\n";
    return 2;
  }

  ofstream dump("configuration.lammpstrj");
  if (!dump) {
    cerr << "could not open configuration.lammpstrj\n";
    return 2;
  }

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
