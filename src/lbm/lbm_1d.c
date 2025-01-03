#include "lbm/lbm_1d.h"

void initialize_1d(lbm_params_1d *lbm) {
    lbm->f = create_2d_array(lbm->NX, lbm->Q);
    lbm->f_new = create_2d_array(lbm->NX, lbm->Q);
    lbm->rho = create_1d_array(lbm->NX);
    lbm->u = create_1d_array(lbm->NX); 

    for (int x = 0; x < lbm->NX; x++) {
        lbm->rho[x] = 2.0;
        lbm->u[x] = 0.;

        for (int k = 0; k < lbm->Q; k++) {
            lbm->f[x][k] = equilibrium_1d(*lbm, k, lbm->rho[x], lbm->u[x]);
        }
    }
}

double equilibrium_1d(lbm_params_1d lbm, int k, double rho, double ux) {
    double cu = lbm.cx[k] * ux;
    double usq = ux * ux;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

void collide_and_stream_1d(lbm_params_1d lbm)
{   
    #pragma omp parallel for
    for (int x = 0; x < lbm.NX; x++) {
        // Compute macroscopic quantities
        lbm.rho[x] = 0.0;
        lbm.u[x] = 0.0;
        for (int k = 0; k < lbm.Q; k++) {
            lbm.rho[x] += lbm.f[x][k];
            lbm.u[x] += lbm.f[x][k] * lbm.cx[k];
        }
        lbm.u[x] /= lbm.rho[x];

        // Collision step
        #pragma omp parallel for
        for (int k = 0; k < lbm.Q; k++) {
            double feq = equilibrium_1d(lbm, k, lbm.rho[x], lbm.u[x]);
            //feq = clamp(feq, .00001, 7.);
            lbm.f[x][k] += (feq - lbm.f[x][k]) / lbm.tau;
            lbm.f[x][k] = fmax(lbm.f[x][k], 0.);
        }
    }

    // Streaming step
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            int xp = x + lbm.cx[k];
            if (xp == -1 || xp >= lbm.NX - 1) {
                lbm.f_new[x][lbm.opp[k]] = lbm.f[x][k];
                //lbm.f[x][k] = 0;
            }
            else lbm.f_new[xp][k] = lbm.f[x][k];
            if (x == lbm.NX - 1) printf_s("B%f B", lbm.f_new[x][k]);
        }       
    }
}

void apply_boundary_conditions(lbm_params_1d lbm) {
    for (int k = 0; k < lbm.Q; k++) {
        lbm.f[0][k] = lbm.f[1][k];
    }
    for (int k = 0; k < lbm.Q; k++) {
        lbm.f[lbm.NX-1][k] = lbm.f[lbm.NX-2][k];
    }
}

void swap_distributions_1d(lbm_params_1d lbm){
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            lbm.f[x][k] = lbm.f_new[x][k];
        }
    }
}

void lbm_1d_step(lbm_params_1d lbm) {
    collide_and_stream_1d(lbm);
    apply_boundary_conditions(lbm);
    swap_distributions_1d(lbm);
}

void lbm_free_1d(lbm_params_1d lbm){
    free_1d_array(lbm.rho);
    free_2d_array(lbm.f, lbm.NX);
    free_2d_array(lbm.f_new, lbm.NX);
    free_1d_array(lbm.u);
}