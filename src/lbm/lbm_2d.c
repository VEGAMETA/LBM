#include "lbm/lbm_2d.h"
#include <string.h>

const int BOXES_NUM = 2;

box BOXES[BOXES_NUM] = {
    {10, 50, 10, 20},
    {30, 40, 50, 65}
}; 

circle CIRCLE = {60, 40, 11};

bool check_inbox(int x, int y){
    for (int i = 0; i < BOXES_NUM; i++)
        if(x >= BOXES[i].x0 && x <= BOXES[i].x1 && y >= BOXES[i].y0 && y <= BOXES[i].y1)
            return 1;
    return 0;
}

bool check_circle(int x, int y, circle _circle){
    if ( (x - _circle.x) * (x - _circle.x) + (y - _circle.y) * (y - _circle.y) <= _circle.r*_circle.r ) return 1;
    return 0;
}

void initialize_2d(lbm_params_2d *lbm) {
    lbm->f = create_3d_array(lbm->NX, lbm->NY,lbm->Q);
    lbm->f_new = create_3d_array(lbm->NX, lbm->NY, lbm->Q);
    lbm->rho = create_2d_array(lbm->NX, lbm->NY);
    lbm->u = create_3d_array(lbm->NX, lbm->NY, lbm->D);
    if (lbm->map) map_initialize_2d(*lbm);
    else base_initialize_2d(*lbm);
}

void map_initialize_2d(lbm_params_2d lbm) {
    for (int x = 0; x < lbm.NX; x++){
        for (int y = 0; y < lbm.NY; y++){
            switch (lbm.map[x][y])
            {
                case ADIABATIC:
                    lbm.rho[x][y] = 0;
                    break;
                case CONVECTIVE_OUT:
                    lbm.rho[x][y] = 2;
                    break;
                default:
                    lbm.rho[x][y] = 1;
            }
            if (lbm.rho[x][y] == 0.) continue;
            lbm.u[x][y][0] = lbm.u[x][y][1] = 0.;
            for (int q = 0; q < lbm.Q; q++)
                lbm.f[x][y][q] = equilibrium_2d(lbm, q, lbm.rho[x][y], lbm.u[x][y][0], lbm.u[x][y][1]);
        }
    }
}

void base_initialize_2d(lbm_params_2d lbm) {
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            lbm.rho[x][y] = 1.0;
            if(check_inbox(x, y))
            {
                lbm.rho[x][y] = 0.;
                continue;
            }
            if (check_circle(x, y, CIRCLE)) lbm.rho[x][y] = 2.;
            lbm.u[x][y][0] = lbm.u[x][y][1] = 0.;
            for (int q = 0; q < lbm.Q; q++) {
                lbm.f[x][y][q] = equilibrium_2d(lbm, q, lbm.rho[x][y], lbm.u[x][y][0], lbm.u[x][y][1]);
            }
        }
    }
}

double equilibrium_2d(lbm_params_2d lbm, int q, double rho, double ux, double uy) {
    return lbm.w[q] * rho;
    double cu = lbm.cx[q] * ux + lbm.cy[q] * uy;
    double usq = ux * ux + uy * uy;
    return lbm.w[q] * rho * (1. + 3. * cu + 4.5 * cu * cu - 1.5 * usq);
}

double get_tau(double rho){
    return (0.25 + 0.25 * rho * rho) * 2;
}

void map_collide_2d(lbm_params_2d lbm) {
    #pragma omp parallel for collapse(2)
    for(int x = 0; x < lbm.NX; x++) {
        for(int y = 0; y < lbm.NY; y++) {
            if (lbm.map[x][y] == ADIABATIC) continue;
            lbm.rho[x][y] = 0.;
            lbm.u[x][y][0] = 0.;
            lbm.u[x][y][1] = 0.;
            for(int q = 0; q < lbm.Q; q++) lbm.rho[x][y] += lbm.f[x][y][q];
            if (lbm.map[x][y] == CONVECTIVE_OUT) lbm.rho[x][y] = 2.;
            if (lbm.map[x][y] == CONVECTIVE_IN) lbm.rho[x][y] = 1.;
            lbm.u[x][y][0] = (lbm.f[x][y][1] - lbm.f[x][y][2] + lbm.f[x][y][5] - lbm.f[x][y][6] + lbm.f[x][y][7] - lbm.f[x][y][8]);
            lbm.u[x][y][1] = (lbm.f[x][y][3] - lbm.f[x][y][4] + lbm.f[x][y][5] - lbm.f[x][y][6] + lbm.f[x][y][8] - lbm.f[x][y][7]);
            lbm.u[x][y][0] /= lbm.rho[x][y];
            lbm.u[x][y][1] /= lbm.rho[x][y];
            #pragma omp parallel for
            for(int q = 0; q < lbm.Q; q++) {
                double feq = equilibrium_2d(lbm, q, lbm.rho[x][y], lbm.u[x][y][0], lbm.u[x][y][1]);
                lbm.f[x][y][q] += (feq - lbm.f[x][y][q]) / get_tau(lbm.rho[x][y]);
            }
        }
    }
}
    
void collide_2d(lbm_params_2d lbm) {
    #pragma omp parallel for collapse(2)
    for(int x = 0; x < lbm.NX; x++) {
        for(int y = 0; y < lbm.NY; y++) {
            if(check_inbox(x, y)) continue;
            lbm.rho[x][y] = 0.;
            lbm.u[x][y][0] = 0.;
            lbm.u[x][y][1] = 0.;
            for(int q = 0; q < lbm.Q; q++) lbm.rho[x][y] += lbm.f[x][y][q];
            if ( (x + 2 - CIRCLE.x) * (x - 2 - CIRCLE.x) + (y + 2 - CIRCLE.y) * (y - 2 - CIRCLE.y) <= CIRCLE.r*CIRCLE.r ) lbm.rho[x][y] = 2.;
            if (y == 0 || x == lbm.NX - 1 ||  x == 0) lbm.rho[x][y] = 1.0;
            for (int i = 0; i < BOXES_NUM; i++){
                if(x >= BOXES[i].x0 - 1 && x <= BOXES[i].x1 + 1 && y >= BOXES[i].y0 - 1 && y <= BOXES[i].y1 + 1)
                    lbm.rho[x][y] = 1.;
            }
            lbm.u[x][y][0] = (lbm.f[x][y][1] - lbm.f[x][y][2] + lbm.f[x][y][5] - lbm.f[x][y][6] + lbm.f[x][y][7] - lbm.f[x][y][8]);
            lbm.u[x][y][1] = (lbm.f[x][y][3] - lbm.f[x][y][4] + lbm.f[x][y][5] - lbm.f[x][y][6] + lbm.f[x][y][8] - lbm.f[x][y][7]);
            lbm.u[x][y][0] /= lbm.rho[x][y];
            lbm.u[x][y][1] /= lbm.rho[x][y];
            #pragma omp parallel for
            for(int q = 0; q < lbm.Q; q++) {
                double feq = equilibrium_2d(lbm, q, lbm.rho[x][y], lbm.u[x][y][0], lbm.u[x][y][1]);
                lbm.f[x][y][q] += (feq - lbm.f[x][y][q]) / get_tau(lbm.rho[x][y]);
            }
        }
    }
}

void map_stream_2d(lbm_params_2d lbm) {
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int q = 0; q < lbm.Q; q++) {
                if (lbm.map[x][y] == ADIABATIC) continue;
                int _x = x + lbm.cx[q];
                int _y = y + lbm.cy[q];
                if (_x == lbm.NX || _y == lbm.NY || _x == -1 || _y == -1)
                    lbm.f_new[x][y][lbm.opp[q]] = lbm.f[x][y][q];
                else if (lbm.map[_x][_y] == ADIABATIC){
                    lbm.f_new[x][y][lbm.opp[q]] = lbm.f[x][y][q];
                    lbm.f_new[_x][_y][lbm.opp[q]] = lbm.f[_x][_y][q];
                }
                else lbm.f_new[_x][_y][q] = lbm.f[x][y][q];
            }
        }
    }
}

void stream_2d(lbm_params_2d lbm) {
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int q = 0; q < lbm.Q; q++) {
                if (check_inbox(x, y)) continue;
                int _x = x + lbm.cx[q];
                int _y = y + lbm.cy[q];
                if (_x == lbm.NX || _y == lbm.NY || _x == -1 || _y == -1)
                    lbm.f_new[x][y][lbm.opp[q]] = lbm.f[x][y][q];
                else if (check_inbox(_x, _y)){
                    lbm.f_new[x][y][lbm.opp[q]] = lbm.f[x][y][q];
                    lbm.f_new[_x][_y][lbm.opp[q]] = lbm.f[_x][_y][q];
                }
                else lbm.f_new[_x][_y][q] = lbm.f[x][y][q];
            }
        }
    }
}

void swap_distributions_2d(lbm_params_2d *lbm){
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < lbm->NX; x++) {
        for (int y = 0; y < lbm->NY; y++) {
            memcpy(lbm->f[x][y], lbm->f_new[x][y], sizeof(double) * lbm->Q);
        }
    }
}

void lbm_2d_step(lbm_params_2d lbm) {
    if (lbm.map){
        map_collide_2d(lbm);
        map_stream_2d(lbm);
    }
    else{
        collide_2d(lbm);
        stream_2d(lbm);
        CIRCLE.y = 40. + 20. * time_sin();
        CIRCLE.x = 45. + 20. * time_cos();
    }
    swap_distributions_2d(&lbm);
}

void lbm_free_2d(lbm_params_2d lbm) {
    free_2d_array(lbm.rho, lbm.NX);
    free_3d_array(lbm.f, lbm.NX, lbm.NY);
    free_3d_array(lbm.f_new, lbm.NX, lbm.NY);
    free_3d_array(lbm.u, lbm.NX, lbm.NY);
}