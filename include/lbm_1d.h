#pragma once
#include "utils.h"

typedef struct lbm_parameters_1d
{
    const int NX;
    int D;
    const int Q;
    const double *w;
    const int *cx;
    const int *opp;
    const double EQ_A;
    const double EQ_B;
    const double EQ_C;
    const double EQ_D;

    double tau;

    double **f;
    double **f_new;
    double *rho;
    double *u;
} lbm_params_1d;

void initialize_1d(lbm_params_1d *lbm);
double equilibrium_1d(lbm_params_1d lbm, int k, double rho, double ux);
void collide_and_stream_1d(lbm_params_1d lbm);
void apply_boundary_conditions(lbm_params_1d lbm);
void swap_distributions_1d(lbm_params_1d lbm);
void lbm_1d_step(lbm_params_1d lbm);
void lbm_free_1d(lbm_params_1d lbm);