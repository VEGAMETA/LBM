#include <stdio.h>
#include <stdlib.h>
#include "lbm/lbm.h"

const int STEPS = 10;

void d3q19_test(int steps) {
    const double w[19] = {
        1.0/3.0,
        1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
    };
    const int cx[19] = { 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1};
    const int cy[19] = { 0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1,  1, -1,  1, -1};
    const int cz[19] = { 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1,  0,  0};
    double tau = 0.5; // Relaxation time
    lbm_params lbm = {50, 50, 50, 3, 19, w, cx, cy, cz, 1.0, 3.0, 4.5, -1.5, tau};

    lbm_3d_loop(lbm, steps);
}

void d2q9_test(int steps) {
    const double w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    const int cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1}; // Velocity directions in x
    const int cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1}; // Velocity directions in y
    double tau = 0.6;
    lbm_params lbm = {50, 50, 0, 3, 9, w, cx, cy, 0, 1.0, 3.0, 4.5, -1.5, tau};

    lbm_2d_loop(lbm, steps);
}

int main() {
    puts("running...");
    d3q19_test(1000);
    return 0;
}
