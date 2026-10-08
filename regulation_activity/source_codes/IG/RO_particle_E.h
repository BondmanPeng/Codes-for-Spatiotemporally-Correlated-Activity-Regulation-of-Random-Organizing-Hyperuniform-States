#ifndef RO_PARTICLE_H
#define RO_PARTICLE_H

#include <time.h>
#include <cmath>
#include <random>
#include <algorithm>

static std::mt19937 rng(time(nullptr));

inline double randCV() {
    return rng() * (1.0 / 4294967296.0);
}

inline unsigned int randNum() {
    return rng();
}

struct Particle {
   int id;
   int type;
   double radius;
   double weight;
    double boxs[2];
    double backup_boxs[2];
    double dE;
    double average_disp;

    Particle* cell_next_pointer;
    Particle* cell_pre_pointer;
    Particle* list_next_pointer;
    Particle* list_pre_pointer;

    bool   active;
    double disp[2];
    double cum_disp[2];
    int    n_overlaps;

    void initialize(int _id, int _type, double _weight, double _radius) {
           id = _id;
           type = _type;
           weight = _weight;
           radius=_radius;
           dE=1.*_weight;
        boxs[0] = randCV();
        boxs[1] = randCV();
        backup_boxs[0] = boxs[0];
        backup_boxs[1] = boxs[1];
        cell_next_pointer  = nullptr;
        cell_pre_pointer   = nullptr;
        list_next_pointer  = nullptr;
        list_pre_pointer   = nullptr;
        active = false;
        disp[0] = 0.0;
        disp[1] = 0.0;
        cum_disp[0] = 0.0;
        cum_disp[1] = 0.0;
        n_overlaps = 0;
    }

    void backup() {
        backup_boxs[0] = boxs[0];
        backup_boxs[1] = boxs[1];
    }

    void reject() {
        boxs[0] = backup_boxs[0];
        boxs[1] = backup_boxs[1];
    }

    void reset_disp() {
        disp[0] = 0.0;
        disp[1] = 0.0;
        n_overlaps = 0;
        active = false;
    }
};

#endif
