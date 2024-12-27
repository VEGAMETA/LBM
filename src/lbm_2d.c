#include "lbm_2d.h"

box BOXES[4] = {
    //{50, 150, 220, 295},
    //{450, 550, 195, 220},
    //{150, 400, 50, 120},
    //{300, 350, 175, 350}
    {0, 250, 300, 400},
    {0, 250, 0, 100}
}; // x0, x1, y0, y1  x0 < x1, y0 < y1

bool check_inbox(int x, int y){
    for (int i = 0; i < 2; i++){
        if(x >= BOXES[i].x0 && x <= BOXES[i].x1 && y >= BOXES[i].y0 && y <= BOXES[i].y1)
            return 1;
    }
    return 0;
}   

void initialize_2d(lbm_params_2d *lbm) {
    lbm->f = create_3d_array(lbm->NX, lbm->NY,lbm->Q);
    lbm->f_new = create_3d_array(lbm->NX, lbm->NY, lbm->Q);
    lbm->rho = create_2d_array(lbm->NX, lbm->NY);;
    lbm->u = create_3d_array(lbm->NX, lbm->NY, lbm->D); 

    for (int x = 0; x < lbm->NX; x++) {
        for (int y = 0; y < lbm->NY; y++) {
            // Initial density
            lbm->rho[x][y] = 4.0;
            // lbm->rho[x][y] += 1.*
            //     (((x-100)*(x-100) +  (y-250)*(y-250)) < 24000) || 
            //     (x > 450 && y > 120 && y < 180); 
            if(check_inbox(x, y))
            {
                lbm->rho[x][y] = 0.;
                continue;
            }
            lbm->u[x][y][0] = 0.;
            lbm->u[x][y][1] = .0;

            for (int k = 0; k < lbm->Q; k++) {
                lbm->f[x][y][k] = equilibrium_2d(*lbm, k, lbm->rho[x][y], lbm->u[x][y][0], lbm->u[x][y][1]);
            }
        }
    }
}

double equilibrium_2d(lbm_params_2d lbm, int k, double rho, double ux, double uy) {
    double cu = lbm.cx[k] * ux + lbm.cy[k] * uy;
    double usq = ux * ux + uy * uy;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

double get_tau(lbm_params_2d lbm, int x, int y) {
    double t = lbm.rho[x][y]/2. - 1.;
    double alpha = 1. + 0.5*(t*t);
    return alpha * 0.5;
}

void collide_and_stream_2d(lbm_params_2d lbm) {
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            if (check_inbox(x, y)) continue;
            lbm.rho[x][y] = 0.;
            for (int k = 0; k < lbm.Q; k++) {
                lbm.rho[x][y] += lbm.f[x][y][k];
                // lbm.u[x][y][0] += lbm.f[x][y][k] * lbm.cx[k];
                // lbm.u[x][y][1] += lbm.f[x][y][k] * lbm.cy[k];
            }
            lbm.u[x][y][0] = lbm.f[x][y][1] - lbm.f[x][y][3] + lbm.f[x][y][5] - lbm.f[x][y][7] + lbm.f[x][y][8] - lbm.f[x][y][6];
            lbm.u[x][y][1] = lbm.f[x][y][2] - lbm.f[x][y][4] + lbm.f[x][y][5] - lbm.f[x][y][7] - lbm.f[x][y][8] + lbm.f[x][y][6];
            // left wall convective
            if (x == 0) lbm.rho[x][y] = 2.;
            if (lbm.rho[x][y] <= 1.){
                lbm.rho[x][y] = 1.;
                lbm.u[x][y][0] = lbm.u[x][y][1] = 0.;
            }
            else{
                lbm.u[x][y][0] /= lbm.rho[x][y];
                lbm.u[x][y][1] /= lbm.rho[x][y];
            }
            

            // Collision step
            #pragma omp parallel for
            for (int k = 0; k < lbm.Q; k++) {
                double feq = equilibrium_2d(lbm, k, lbm.rho[x][y], lbm.u[x][y][0], lbm.u[x][y][1]);
                //feq = clamp(feq, .00001, 7.);
                lbm.f[x][y][k] += (feq - lbm.f[x][y][k]) / get_tau(lbm, x, y); ///lbm.tau;
                lbm.f[x][y][k] = fmax(lbm.f[x][y][k], 0.);
            }
        }
    }

    // Streaming step
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int k = 0; k < lbm.Q; k++) {
                if (check_inbox(x, y)) continue;
                // int xp = (x + lbm.cx[k] + lbm.NX) % lbm.NX;
                // int yp = (y + lbm.cy[k] + lbm.NY) % lbm.NY;
                // if (xp == 0 || xp == lbm.NX - 1 || yp == 0 || yp == lbm.NY - 1)
                //     lbm.f_new[xp][yp][lbm.opp[k]] = 3.; // convective
                int xp = x + lbm.cx[k];
                int yp = y + lbm.cy[k];
                if (yp == - 1 || yp == lbm.NY || xp == - 1 || xp == lbm.NX)
                    lbm.f_new[x][y][lbm.opp[k]] = lbm.f[x][y][k];//adiobatic
                else if (check_inbox(xp, yp)) {
                    lbm.f_new[xp][yp][lbm.opp[k]] = lbm.f[xp][yp][k];
                    lbm.f_new[x][y][lbm.opp[k]] = lbm.f[x][y][k];
                }
                else lbm.f_new[xp][yp][k] = lbm.f[x][y][k];
            }
        }
    }
}


void boundary(lbm_params_2d *lbm, int x0, int x1, int y0, int y1) {
    for (int x = x0 + 1; x < x1; x++){
        for (int q = 0; q < lbm->Q; q++){
            lbm->f_new[x][y0][lbm->opp[q]] = lbm->f[x][y0][q];
            lbm->f_new[x][y1][lbm->opp[q]] = lbm->f[x][y1][q];
        }
    }
    for (int y = y0 + 1; y < y1; y++){
        for (int q = 0; q < lbm->Q; q++){
            lbm->f_new[x0][y][lbm->opp_x[q]] = lbm->f[x0][y][q];
            lbm->f_new[x1][y][lbm->opp_x[q]] = lbm->f[x1][y][q];
        }
    }
    
    for (int q = 0; q < lbm->Q; q++){
        lbm->f_new[x0][y0][lbm->opp[q]] = lbm->f[x0][y0][q];
        lbm->f_new[x1][y1][lbm->opp[q]] = lbm->f[x1][y1][q];
        lbm->f_new[x0][y1][lbm->opp[q]] = lbm->f[x0][y1][q];
        lbm->f_new[x1][y0][lbm->opp[q]] = lbm->f[x1][y0][q];
    }
}


void apply_boundary_conditions_2d(lbm_params_2d *lbm) {
    // boundary(lbm, 50, 150, 220, 295);
    // boundary(lbm, 450, 550, 195, 220);
    // boundary(lbm, 150, 400, 50, 120);
    // boundary(lbm, 300, 350, 175, 350);
}

void swap_distributions_2d(lbm_params_2d lbm){
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int k = 0; k < lbm.Q; k++) {
                lbm.f[x][y][k] = lbm.f_new[x][y][k];
            }
        }
    }
}

void lbm_2d_step(lbm_params_2d lbm) {
    collide_and_stream_2d(lbm);
    apply_boundary_conditions_2d(&lbm);
    swap_distributions_2d(lbm);
}

void lbm_free_2d(lbm_params_2d lbm) {
    free_2d_array(lbm.rho, lbm.NX);
    free_3d_array(lbm.f, lbm.NX, lbm.NY);
    free_3d_array(lbm.f_new, lbm.NX, lbm.NY);
    free_3d_array(lbm.u, lbm.NX, lbm.NY);
}
