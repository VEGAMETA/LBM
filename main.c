#include <stdio.h>
#include <stdlib.h>
#include "lbm.h"
#include "visualize_raylib.h"

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
    double tau = 0.8; // Relaxation time
    lbm_params_3d lbm = {50, 50, 50, 3, 19, w, cx, cy, cz, 1.0, 3.0, 4.5, -1.5, tau};

}

void d2q9_test() {
    const double w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    const int cx[9] = {0, 1, 0, -1,  0, 1, -1, -1,  1}; // Velocity directions in x
    const int cy[9] = {0, 0, 1,  0, -1, 1,  1, -1, -1}; // Velocity directions in y
    // 6 2 5  ^y
    //  \|/   |   x
    // 3-0-1   --->
    //  /|\
    // 7 4 8
    int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    int opp_x[9] = {0, 3, 4, 1, 2, 6, 5, 8, 7};
    int opp_y[9] = {0, 3, 4, 1, 2, 8, 7, 6, 5};
    double tau = .65;
    lbm_params_2d lbm = {450, 400, 3, 9, w, cx, cy, opp, opp_x, opp_y, 1.0, 3.0, 4.5, -1.5, tau};
    initialize_2d(&lbm);
    visualize_2d(lbm);
    lbm_free_2d(lbm);
}

void d1q2_test() {
    const double w[2] = {1.0/2.0, 1.0/2.0};
    const int cx[2] = {-1, 1};
    const int opp[2] = {1, -1};
    double tau = .75;
    lbm_params_1d lbm = {600, 3, 2, w, cx, opp, 1.0, 3.0, 4.5, -1.5, tau};
    initialize_1d(&lbm);
    visualize_1d(lbm);
    lbm_free_1d(lbm);
}

int main(int argc, char* argv[]) {
    int steps = atoi(argv[1]); 
    //printf("running %d steps", steps); 
    d2q9_test();
    //d1q2_test();
    return 0;
}
