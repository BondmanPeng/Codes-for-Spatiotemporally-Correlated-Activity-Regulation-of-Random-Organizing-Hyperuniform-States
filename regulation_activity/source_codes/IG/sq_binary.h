#pragma once

#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

struct Atom { int id; int type; double x, y, vx, vy, radius; };
struct Box {
  double xlo, xhi, ylo, yhi;
  double lx() const { return xhi - xlo; }
  double ly() const { return yhi - ylo; }
};

inline std::vector<Atom> read_current_frame(std::ifstream &file, Box &box) {
  std::string line; std::vector<Atom> atoms;
  if (!std::getline(file, line)) return atoms;
  if (!std::getline(file, line)) return atoms;
  if (!std::getline(file, line)) return atoms;
  if (!std::getline(file, line)) return atoms;
  int num_atoms = 0; { std::stringstream ss(line); ss >> num_atoms; }
  if (num_atoms <= 0) return atoms;
  atoms.reserve(num_atoms);
  if (!std::getline(file, line)) return atoms;
  std::stringstream ss; double wlo, whi;
  if (!std::getline(file, line)) return atoms; ss.clear(); ss.str(line); ss >> box.xlo >> box.xhi;
  if (!std::getline(file, line)) return atoms; ss.clear(); ss.str(line); ss >> box.ylo >> box.yhi;
  if (!std::getline(file, line)) return atoms; ss.clear(); ss.str(line); ss >> wlo >> whi;
  if (!std::getline(file, line)) return atoms;
  for (int i = 0; i < num_atoms; ++i) {
    if (!std::getline(file, line)) { atoms.clear(); return atoms; }
    std::stringstream ads(line); Atom a; double wz, wvz;
    ads >> a.id >> a.type >> a.x >> a.y >> wz >> a.vx >> a.vy >> wvz >> a.radius;
    atoms.push_back(a);
  }
  return atoms;
}

struct SqAccumulator {
  int max_mode = 0; double dk = 0.0; int frames = 0;
  std::vector<double> q, S_AA, S_BB, S_AB, S_nn;
  std::vector<int> shell_counts;
  void init(int mm, double dk_) {
    max_mode = mm; dk = dk_; frames = 0;
    q.assign(mm, 0.0); S_AA.assign(mm, 0.0); S_BB.assign(mm, 0.0);
    S_AB.assign(mm, 0.0); S_nn.assign(mm, 0.0); shell_counts.assign(mm, 0);
  }
};

struct ShellPartials {
  std::vector<double> q, aa, bb, ab, nn; std::vector<long> cnt;
  void init(int n) {
    q.assign(n, 0.0); aa.assign(n, 0.0); bb.assign(n, 0.0);
    ab.assign(n, 0.0); nn.assign(n, 0.0); cnt.assign(n, 0);
  }
};

inline void accumulate_frame_sq(SqAccumulator &acc,
                                const std::vector<Atom> &atoms, const Box &box) {
  const int max_mode = acc.max_mode, N = (int)atoms.size();
  const double dk = acc.dk;
  if (N == 0 || max_mode <= 0) return;

  int NA = 0, NB = 0;
  std::vector<unsigned char> isA(N);
  for (int i = 0; i < N; ++i) { isA[i] = (atoms[i].type == 1); isA[i] ? ++NA : ++NB; }
  const double inv_NA = 1.0 / (NA ? NA : 1), inv_NB = 1.0 / (NB ? NB : 1);
  const double inv_NAB = 1.0 / std::sqrt((double)(NA ? NA : 1) * (NB ? NB : 1));
  const double inv_N = 1.0 / N;

  const int M1 = max_mode + 1;
  const std::size_t str = (std::size_t)N;
  std::vector<double> CX((std::size_t)M1 * str), SX((std::size_t)M1 * str),
      CY((std::size_t)M1 * str), SY((std::size_t)M1 * str);
  {
    std::vector<double> cdx(N), sdx(N), cdy(N), sdy(N);
    for (int i = 0; i < N; ++i) {
      CX[i] = 1; SX[i] = 0; CY[i] = 1; SY[i] = 0;
      cdx[i] = std::cos(dk * atoms[i].x); sdx[i] = std::sin(dk * atoms[i].x);
      cdy[i] = std::cos(dk * atoms[i].y); sdy[i] = std::sin(dk * atoms[i].y);
    }
    for (int m = 1; m < M1; ++m) {
      const std::size_t c = (std::size_t)m * str, p = (std::size_t)(m - 1) * str;
      for (int i = 0; i < N; ++i) {
        CX[c + i] = CX[p + i] * cdx[i] - SX[p + i] * sdx[i];
        SX[c + i] = SX[p + i] * cdx[i] + CX[p + i] * sdx[i];
        CY[c + i] = CY[p + i] * cdy[i] - SY[p + i] * sdy[i];
        SY[c + i] = SY[p + i] * cdy[i] + CY[p + i] * sdy[i];
      }
    }
  }

#ifdef _OPENMP
  const int nth = omp_get_max_threads();
#else
  const int nth = 1;
#endif
  std::vector<ShellPartials> tls(nth);
  for (auto &t : tls) t.init(max_mode);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    ShellPartials &sa = tls[tid];
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1) nowait
#endif
    for (int m = 0; m <= max_mode; ++m) {
      const std::size_t mx = (std::size_t)m * str;
      const int n_lo = (m == 0) ? 1 : -max_mode;
      for (int n = n_lo; n <= max_mode; ++n) {
        const int mn2 = m * m + n * n;
        const int slot = (int)std::floor(std::sqrt((double)mn2)) - 1;
        if (slot < 0 || slot >= max_mode) continue;
        const int an = n < 0 ? -n : n;
        const double nsg = n < 0 ? -1.0 : 1.0;
        const std::size_t ny = (std::size_t)an * str;

        double Ar = 0, Ai = 0, Br = 0, Bi = 0;
        for (int i = 0; i < N; ++i) {
          const double cx = CX[mx + i], sx = SX[mx + i];
          const double cy = CY[ny + i], sy = nsg * SY[ny + i];
          const double c = cx * cy - sx * sy;
          const double s = sx * cy + cx * sy;
          if (isA[i]) { Ar += c; Ai += s; } else { Br += c; Bi += s; }
        }
        const double nr = Ar + Br, ni = Ai + Bi;
        const double w = 2.0;
        sa.q[slot]  += w * std::sqrt((double)mn2) * dk;
        sa.aa[slot] += w * (Ar * Ar + Ai * Ai) * inv_NA;
        sa.bb[slot] += w * (Br * Br + Bi * Bi) * inv_NB;
        sa.ab[slot] += w * (Ar * Br + Ai * Bi) * inv_NAB;
        sa.nn[slot] += w * (nr * nr + ni * ni) * inv_N;
        sa.cnt[slot] += 2;
      }
    }
  }

  std::vector<double> q_s(max_mode, 0), aa(max_mode, 0), bb(max_mode, 0),
      ab(max_mode, 0), nn(max_mode, 0);
  std::vector<long> cnt(max_mode, 0);
  for (int t = 0; t < nth; ++t) {
    const ShellPartials &sa = tls[t];
    for (int s = 0; s < max_mode; ++s) {
      q_s[s] += sa.q[s]; aa[s] += sa.aa[s]; bb[s] += sa.bb[s];
      ab[s] += sa.ab[s]; nn[s] += sa.nn[s]; cnt[s] += sa.cnt[s];
    }
  }
  const double prev = acc.frames, curr = prev + 1.0;
  for (int s = 0; s < max_mode; ++s) {
    if (!cnt[s]) continue;
    const double ic = 1.0 / (double)cnt[s];
    acc.q[s]    = (acc.q[s]    * prev + q_s[s] * ic) / curr;
    acc.S_AA[s] = (acc.S_AA[s] * prev + aa[s]  * ic) / curr;
    acc.S_BB[s] = (acc.S_BB[s] * prev + bb[s]  * ic) / curr;
    acc.S_AB[s] = (acc.S_AB[s] * prev + ab[s]  * ic) / curr;
    acc.S_nn[s] = (acc.S_nn[s] * prev + nn[s]  * ic) / curr;
    acc.shell_counts[s] = (int)cnt[s];
  }
  acc.frames++;
}
