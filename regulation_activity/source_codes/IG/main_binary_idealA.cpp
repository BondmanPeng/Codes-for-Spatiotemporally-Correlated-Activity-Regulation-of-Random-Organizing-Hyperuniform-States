#include "randorg_binary_idealA.h"

using namespace std;

int main(int argc, char *argv[]) {
  ofstream dump("configuration.lammpstrj");

  double L = (argc > 9) ? stod(argv[9]) : 300.0;
  double radius_A = 0.5, radius_B = 0.5;
  double rho_A = (argc > 1) ? stod(argv[1]) : 0.60;
  double rho_B = (argc > 2) ? stod(argv[2]) : 0.60;
  double h = (argc > 3) ? stod(argv[3]) : 1.00;
  double h0 = (argc > 4) ? stod(argv[4]) : 2.0;
  double T = (argc > 5) ? stod(argv[5]) : 0.10;
  double R_A = (argc > 6) ? stod(argv[6]) : 5.00;
  int steps = (argc > 7) ? stoi(argv[7]) : 20000;
  int n_snap = (argc > 8) ? stoi(argv[8]) : 200;

  if (L <= 0.0) {
    cerr << "Require L > 0\n";
    return 1;
  }

  int dump_start = steps / 2;
  int interval = max(1, (steps - dump_start) / n_snap);

  RandOrgBinary sys;
  sys.initialize_system(L, radius_A, radius_B, rho_A, rho_B, h, h0, R_A, T);
  printf("  mean A number density = %.4f\n", sys.mean_A_number_density);
  printf("  dumping %d snapshots every %d steps from t=%d\n", n_snap, interval,
         dump_start);

  sys.run(steps, dump, interval, dump_start,             true);
  dump.close();
  sys.cleanup();

  printf("\nSaved configuration.lammpstrj. Structure factor: run "
         "compute_sq_species\n");
  return 0;
}
