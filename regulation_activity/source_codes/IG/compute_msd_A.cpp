#include "randorg_binary_idealA.h"

#include <algorithm>
#include <iomanip>

using namespace std;

double mean_squared_displacement_A(const RandOrgBinary &sys) {
  if (sys.N_A == 0)
    return 0.0;

  double sum = 0.0;
  for (Particle *p = sys.listA.headPtr; p; p = p->list_next_pointer) {
    double dx = p->cum_disp[0] * sys.L;
    double dy = p->cum_disp[1] * sys.L;
    sum += dx * dx + dy * dy;
  }
  return sum / sys.N_A;
}

int main(int argc, char *argv[]) {
  double L = 300.0;
  double radius_A = 0.5, radius_B = 0.5;
  double rho_A = (argc > 1) ? stod(argv[1]) : 0.60;
  double rho_B = (argc > 2) ? stod(argv[2]) : 0.60;
  double h = (argc > 3) ? stod(argv[3]) : 1.00;
  double h0 = (argc > 4) ? stod(argv[4]) : 1.00;
  double T = (argc > 5) ? stod(argv[5]) : 0.10;
  double R_A = (argc > 6) ? stod(argv[6]) : 5.00;
  int steps = (argc > 7) ? stoi(argv[7]) : 20000;
  int sample_interval = (argc > 8) ? max(1, stoi(argv[8])) : 100;

  RandOrgBinary sys;
  sys.initialize_system(L, radius_A, radius_B, rho_A, rho_B, h, h0, R_A, T);

  ofstream output("msd_A.csv");
  if (!output) {
    cerr << "cannot open msd_A.csv for writing\n";
    sys.cleanup();
    return 1;
  }

  output << "step,msd_A\n";
  output << "0,0\n";
  output << setprecision(16);

  for (int step = 1; step <= steps; step++) {
    sys.step();
    if (step % sample_interval == 0 || step == steps)
      output << step << ',' << mean_squared_displacement_A(sys) << '\n';
    if (step % 500 == 0)
      printf("  t=%6d  MSD_A=%.8g\n", step,
             mean_squared_displacement_A(sys));
  }

  output.close();
  sys.cleanup();
  printf("\nSaved msd_A.csv\n");
  return 0;
}
