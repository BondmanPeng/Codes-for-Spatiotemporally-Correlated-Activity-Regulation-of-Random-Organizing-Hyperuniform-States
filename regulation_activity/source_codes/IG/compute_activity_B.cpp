#include "randorg_binary_idealA.h"
#include <iomanip>

int main(int argc, char *argv[]) {
  const double rho_A = argc > 1 ? stod(argv[1]) : 0.6;
  const double rho_B = argc > 2 ? stod(argv[2]) : 0.6;
  const double h0 = argc > 3 ? stod(argv[3]) : 2.0;
  const int steps = argc > 4 ? stoi(argv[4]) : 100000;
  const unsigned seed = argc > 5 ? stoul(argv[5]) : 1;
  const double L = argc > 6 ? stod(argv[6]) : 300.0;
  if (!isfinite(rho_A) || !isfinite(rho_B) || !isfinite(h0) ||
      !isfinite(L) || rho_A < 0 || rho_B <= 0 || h0 < 0 ||
      steps < 2 || L < 6 || rho_B * L * L < 1) {
    cerr << "Require finite densities/strength, rho_A>=0, N_B>=1, h0>=0, steps>=2, L>=6\n";
    return 1;
  }
  rng.seed(seed);
  RandOrgBinary sys;
  sys.initialize_system(L, 0.5, 0.5, rho_A, rho_B, 0.0, h0, 3.0, 10.0);
  ofstream activity("activity_B.csv"), summary("activity_summary.csv");
  if (!activity || !summary) {
    cerr << "Cannot open activity output\n";
    sys.cleanup();
    return 1;
  }
  activity << "step,active_count_B,active_fraction_B\n" << setprecision(16);
  double tail_sum = 0.0;
  int last_step = 0, absorbing_step = -1;
  string status = "active_at_limit";
  if (h0 == 0.0) {

    status = "zero_kick";
    activity << "0,0,0\n";
  } else {
    for (int n = 1; n <= steps; ++n) {
      int active = sys.step();
      last_step = n;
      if (n > steps / 2)
        tail_sum += sys.active_frac_B;
      if (n == 1 || n % 100 == 0 || n == steps || active == 0) {
        activity << n << ',' << active << ',' << sys.active_frac_B << '\n';
        activity.flush();
      }
      if (active == 0) {
        absorbing_step = n;
        status = "absorbed";
        break;
      }
    }
  }

  summary << "density_convention,rho_A,rho_B,N_A,N_B,L,h,h0,T,R_A,seed,steps_requested,"
             "steps_run,status,absorbing_step,tail_mean_active_B,final_active_B\n";
  summary << setprecision(16) << "number," << rho_A << ',' << rho_B << ','
          << sys.N_A << ',' << sys.N_B << ',' << L << ",0," << h0 << ",10,3,"
          << seed << ',' << steps << ',' << last_step << ',' << status << ','
          << absorbing_step << ',' << tail_sum / (steps - steps / 2) << ','
          << sys.active_frac_B << '\n';
  cout << "status=" << status << " h0=" << h0 << " seed=" << seed
       << " tail_mean_active_B=" << tail_sum / (steps - steps / 2) << '\n';
  sys.cleanup();
}
