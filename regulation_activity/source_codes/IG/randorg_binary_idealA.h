#ifndef RANDORG_BINARY_IDEALA_H
#define RANDORG_BINARY_IDEALA_H

#include "RO_particle_E.h"
#include "array_computation.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

inline double randGauss() {
  static thread_local std::normal_distribution<double> nd(0.0, 1.0);
  return nd(rng);
}

class Cell_List {
public:
  Particle *headPtr = nullptr;
  void insert(Particle *p) {
    if (headPtr) {
      p->cell_next_pointer = headPtr;
      headPtr->cell_pre_pointer = p;
    } else
      p->cell_next_pointer = nullptr;
    p->cell_pre_pointer = nullptr;
    headPtr = p;
  }
  void remove(Particle *p) {
    if (p != headPtr)
      (p->cell_pre_pointer)->cell_next_pointer = p->cell_next_pointer;
    else
      headPtr = p->cell_next_pointer;
    if (p->cell_next_pointer)
      (p->cell_next_pointer)->cell_pre_pointer = p->cell_pre_pointer;
    p->cell_next_pointer = nullptr;
    p->cell_pre_pointer = nullptr;
  }
};

class Particle_List {
public:
  Particle *headPtr = nullptr;
  void insert(Particle *p) {
    if (headPtr) {
      p->list_next_pointer = headPtr;
      headPtr->list_pre_pointer = p;
    } else
      p->list_next_pointer = nullptr;
    p->list_pre_pointer = nullptr;
    headPtr = p;
  }
};

class RandOrgBinary {
public:

  double L;
  double radius_A, radius_B;
  double rho_A, rho_B;
  double h;
  double h0;
  double R_A;
  double T;
  int N_A = 0, N_B = 0;

  double mean_A_number_density = 0.0;
  double thermal_step = 0.0;

  int cell_num;
  double cell_length;
  Particle_List listA, listB;
  vector<vector<Cell_List>> cellsA, cellsB;

  vector<pair<int, int>> disk_stencil;

  int next_id = 0;
  int step_count = 0;
  double active_frac_A = 0.0, active_frac_B = 0.0;

  void initialize_system(double L_, double radius_A_, double radius_B_,
                         double rho_A_, double rho_B_, double h_, double h0_,
                         double R_A_, double T_) {
    L = L_;
    radius_A = radius_A_;
    radius_B = radius_B_;
    rho_A = rho_A_;
    rho_B = rho_B_;
    h = h_;
    h0 = h0_;
    R_A = R_A_;
    T = T_;

    N_A = (int)(rho_A * L * L);
    N_B = (int)(rho_B * L * L);
    mean_A_number_density = (double)N_A / (L * L);

    thermal_step = sqrt(T) * (2.0 * radius_A) / L;

    double max_diam = 2.0 * max(radius_A, radius_B);
    int cal = (int)floor(L / (max_diam + 0.1));
    cell_num = max(cal, 3);
    cell_length = 1.0 / (double)cell_num;

    cellsA.assign(cell_num, vector<Cell_List>(cell_num));
    cellsB.assign(cell_num, vector<Cell_List>(cell_num));

    double cell_phys = cell_length * L;
    int ring = (int)ceil(R_A / cell_phys);
    int min_offset = -min(ring, cell_num / 2);
    int max_offset = min(ring, (cell_num - 1) / 2);
    disk_stencil.clear();
    for (int dx = min_offset; dx <= max_offset; dx++)
      for (int dy = min_offset; dy <= max_offset; dy++)
        disk_stencil.push_back({dx, dy});

    for (int i = 0; i < N_A; i++)
      spawn(1, randCV(), randCV());
    for (int i = 0; i < N_B; i++)
      spawn(2, randCV(), randCV());

    printf("  Binary ideal-gas A (number densities): N_A=%d (rho_A=%.3f), N_B=%d (rho_B=%.3f), "
           "L=%.1f, cell_num=%d, R_A=%.2f, mean_n_A=%.4f, h=%.3f, "
           "h0=%.3f, T=%.3f (sigma_step=%.4f box)\n",
           N_A, rho_A, N_B, rho_B, L, cell_num, R_A,
           mean_A_number_density, h, h0, T, thermal_step);
  }

  Particle *spawn(int type, double bx, double by) {
    Particle *p = new Particle;
    p->initialize(next_id++, type, 1.0, (type == 1) ? radius_A : radius_B);
    p->boxs[0] = bx;
    p->boxs[1] = by;
    if (type == 1) {
      listA.insert(p);
      put(p, cellsA);
    } else {
      listB.insert(p);
      put(p, cellsB);
    }
    return p;
  }

  void put(Particle *p, vector<vector<Cell_List>> &grid) {
    grid[cell_of(p->boxs[0])][cell_of(p->boxs[1])].insert(p);
  }
  void extract(Particle *p, vector<vector<Cell_List>> &grid) {
    grid[cell_of(p->boxs[0])][cell_of(p->boxs[1])].remove(p);
  }
  int cell_of(double x) const {
    int c = (int)floor(x / cell_length);
    if (c < 0)
      c = 0;
    else if (c >= cell_num)
      c = cell_num - 1;
    return c;
  }

  vector<pair<Particle *, Particle *>>
  find_pairs(Particle *head, vector<vector<Cell_List>> &grid) {
    vector<pair<Particle *, Particle *>> pairs;
    static const int off[8][2] = {{1, -1}, {1, 0},   {1, 1},  {0, -1},
                                  {0, 1},  {-1, -1}, {-1, 0}, {-1, 1}};
    for (Particle *ptc = head; ptc; ptc = ptc->list_next_pointer) {
      int cx = cell_of(ptc->boxs[0]), cy = cell_of(ptc->boxs[1]);
      double dr[2];

      for (Particle *it = grid[cx][cy].headPtr; it;
           it = it->cell_next_pointer) {
        if (it->id <= ptc->id)
          continue;
        Substract(dr, it->boxs, ptc->boxs);
        if (Norm(dr) * L < ptc->radius + it->radius)
          pairs.push_back({ptc, it});
      }

      for (auto &o : off) {
        int nx = cx + o[0], ny = cy + o[1], adj[2] = {0, 0};
        if (nx < 0)
          adj[0] = -1;
        else if (nx >= cell_num)
          adj[0] = 1;
        if (ny < 0)
          adj[1] = -1;
        else if (ny >= cell_num)
          adj[1] = 1;
        nx -= cell_num * adj[0];
        ny -= cell_num * adj[1];
        for (Particle *it = grid[nx][ny].headPtr; it;
             it = it->cell_next_pointer) {
          if (it->id <= ptc->id)
            continue;
          double img[2];
          periodic_boundary(img, it->boxs, adj);
          Substract(dr, img, ptc->boxs);
          if (Norm(dr) * L < ptc->radius + it->radius)
            pairs.push_back({ptc, it});
        }
      }
    }
    return pairs;
  }

  double A_number_density(const double mid[2]) const {
    int cx = cell_of(mid[0]), cy = cell_of(mid[1]);
    int count = 0;
    double radius_box = R_A / L;
    double radius_box_sq = radius_box * radius_box;
    for (const auto &[dx, dy] : disk_stencil) {
      int nx = cx + dx;
      if (nx < 0)
        nx += cell_num;
      else if (nx >= cell_num)
        nx -= cell_num;
      int ny = cy + dy;
      if (ny < 0)
        ny += cell_num;
      else if (ny >= cell_num)
        ny -= cell_num;
      for (Particle *p = cellsA[nx][ny].headPtr; p;
           p = p->cell_next_pointer) {
        double drx = p->boxs[0] - mid[0];
        double dry = p->boxs[1] - mid[1];
        drx -= round(drx);
        dry -= round(dry);
        if (drx * drx + dry * dry <= radius_box_sq)
          count++;
      }
    }
    return count / (M_PI * R_A * R_A);
  }

  static void pair_geometry(Particle *pi, Particle *pj, double axis[2],
                            double mid[2], double &dist_box) {
    double dr[2] = {pj->boxs[0] - pi->boxs[0], pj->boxs[1] - pi->boxs[1]};
    dr[0] -= round(dr[0]);
    dr[1] -= round(dr[1]);
    dist_box = Norm(dr);
    if (dist_box < 1e-12) {
      double a = randCV() * 2.0 * M_PI;
      axis[0] = cos(a);
      axis[1] = sin(a);
    } else {
      axis[0] = dr[0] / dist_box;
      axis[1] = dr[1] / dist_box;
    }
    mid[0] = pi->boxs[0] + 0.5 * dr[0];
    mid[1] = pi->boxs[1] + 0.5 * dr[1];
    mid[0] -= floor(mid[0]);
    mid[1] -= floor(mid[1]);
  }

  int step() {
    auto pairsB = find_pairs(listB.headPtr, cellsB);

    vector<Particle *> activeA;
    for (Particle *p = listA.headPtr; p; p = p->list_next_pointer) {
      p->disp[0] = thermal_step * randGauss();
      p->disp[1] = thermal_step * randGauss();
      activeA.push_back(p);
    }

    vector<Particle *> activeB;
    double sigma_B = 2.0 * radius_B;
    for (auto &[pi, pj] : pairsB) {
      double axis[2], mid[2], dist_box;
      pair_geometry(pi, pj, axis, mid, dist_box);
      double strength = h0;
      if (h != 0.0)
        strength += h * (A_number_density(mid) - mean_A_number_density);
      if (strength == 0.0)
        continue;
      double db = (randCV() * strength * sigma_B) / L;
      pi->disp[0] -= db * axis[0];
      pi->disp[1] -= db * axis[1];
      pj->disp[0] += db * axis[0];
      pj->disp[1] += db * axis[1];
      if (!pi->active) {
        pi->active = true;
        activeB.push_back(pi);
      }
      if (!pj->active) {
        pj->active = true;
        activeB.push_back(pj);
      }
    }

    apply(activeA, cellsA);
    apply(activeB, cellsB);

    active_frac_A = 1.0;
    active_frac_B = N_B ? (double)activeB.size() / N_B : 0.0;
    step_count++;

    return (int)activeB.size();
  }

  void apply(vector<Particle *> &active, vector<vector<Cell_List>> &grid) {
    for (Particle *p : active) {
      extract(p, grid);
      p->boxs[0] += p->disp[0];
      p->boxs[1] += p->disp[1];
      p->cum_disp[0] += p->disp[0];
      p->cum_disp[1] += p->disp[1];
      p->boxs[0] -= floor(p->boxs[0]);
      p->boxs[1] -= floor(p->boxs[1]);
      put(p, grid);
      p->reset_disp();
    }
  }

  int run(int max_steps, ofstream &dump, int interval, int dump_start = 0,
          bool verbose = true) {
    for (int t = 0; t < max_steps; t++) {
      int na = step();
      if (verbose && t % 500 == 0)
        printf("  t=%6d  active_A=%.4f  active_B=%.4f\n", t, active_frac_A,
               active_frac_B);
      if (t >= dump_start && (t - dump_start) % interval == 0)
        record_configuration(dump, t);
      if (na == 0) {
        if (verbose)
          printf("  B absorbed at step %d\n", t);
        return t;
      }
    }
    return max_steps;
  }

  void record_configuration(ofstream &f, int timestep) {
    f << "ITEM: TIMESTEP\n" << timestep << "\n";
    f << "ITEM: NUMBER OF ATOMS\n" << (N_A + N_B) << "\n";
    f << "ITEM: BOX BOUNDS pp pp pp\n";
    f << "0.0 " << L << "\n" << "0.0 " << L << "\n" << "0.0 0.0\n";
    f << "ITEM: ATOMS id type x y z vx vy vz rad\n";
    char buf[96];
    for (Particle *head : {listA.headPtr, listB.headPtr})
      for (Particle *p = head; p; p = p->list_next_pointer) {
        int n =
            snprintf(buf, sizeof buf, "%d %d %.5f %.5f 0 0 0 0 %.5f\n", p->id,
                     p->type, p->boxs[0] * L, p->boxs[1] * L, p->radius);
        f.write(buf, n);
      }
  }

  void cleanup() {
    for (auto *head : {listA.headPtr, listB.headPtr})
      for (Particle *p = head; p;) {
        Particle *n = p->list_next_pointer;
        delete p;
        p = n;
      }
    listA.headPtr = listB.headPtr = nullptr;
  }
};

#endif
