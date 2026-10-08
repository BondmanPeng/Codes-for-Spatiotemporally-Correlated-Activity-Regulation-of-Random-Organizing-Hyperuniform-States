#include "sq_binary.h"
#include <iostream>
#include <iomanip>
using namespace std;

int main(int argc, char** argv) {
    string fname  = (argc > 1) ? argv[1] : "configuration.lammpstrj";
    double q_max  = (argc > 2) ? stod(argv[2]) : 1.0;
    long   min_ts = (argc > 3) ? stol(argv[3]) : -1;

    ifstream file(fname);
    if (!file.is_open()) { cerr << "cannot open " << fname << "\n"; return 1; }

    vector<streampos> off; vector<long> ts; string line;
    while (true) {
        streampos pos = file.tellg();
        if (!getline(file, line)) break;
        if (line.find("ITEM: TIMESTEP") != string::npos) {
            off.push_back(pos);
            string tl; getline(file, tl); ts.push_back(stol(tl));
        }
    }
    if (off.empty()) { cerr << "no frames\n"; return 1; }
    if (min_ts < 0) min_ts = ts.back() / 2;

    vector<int> sel;
    for (size_t i = 0; i < ts.size(); ++i) if (ts[i] >= min_ts) sel.push_back((int)i);
    cout << "frames total=" << off.size() << ", selected (t>=" << min_ts
         << ")=" << sel.size() << "\n";
    if (sel.empty()) { cerr << "no frames after min_timestep\n"; return 1; }

    SqAccumulator acc; Box box{}; bool ready = false;
    int done = 0;
    for (int idx : sel) {
        file.clear(); file.seekg(off[idx]);
        auto atoms = read_current_frame(file, box);
        if (atoms.empty()) { cerr << "\nframe " << idx << " truncated; skipped\n"; continue; }
        if (!ready) {
            double dk = 2.0 * M_PI / box.lx();
            int mm = (int)floor(q_max / dk);
            if (mm < 1) { cerr << "q_max too small (dk=" << dk << ")\n"; return 1; }
            acc.init(mm, dk);
            cout << "L=" << box.lx() << ", dk=" << dk << ", shells=" << mm
                 << " (q up to " << mm * dk << ")\n";
            ready = true;
        }
        accumulate_frame_sq(acc, atoms, box);
        cout << "\rprocessed " << ++done << "/" << sel.size() << flush;
    }
    cout << "\n";

    ofstream out("sq_species.csv");
    out << "# frames " << acc.frames << " q_max " << q_max << " dk " << acc.dk << "\n";
    out << "q,S_AA,S_BB,S_AB,S_nn,k_count\n";
    out << scientific << setprecision(8);
    for (int s = 0; s < acc.max_mode; ++s) {
        if (!acc.shell_counts[s]) continue;
        out << acc.q[s] << "," << acc.S_AA[s] << "," << acc.S_BB[s] << ","
            << acc.S_AB[s] << "," << acc.S_nn[s] << "," << acc.shell_counts[s] << "\n";
    }
    cout << "wrote sq_species.csv (" << acc.frames << " frames)\n";
    return 0;
}
