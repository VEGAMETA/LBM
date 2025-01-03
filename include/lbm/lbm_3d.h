#pragma once
#include "lbm/utils.h"

typedef struct lbm_parameters_3d
{
    const int NX;
    const int NY;
    const int NZ;
    int D;
    const int Q;

    const double *w;

    const int *cx;
    const int *cy;
    const int *cz;

    const double EQ_A;
    const double EQ_B;
    const double EQ_C;
    const double EQ_D;

    double tau;

    double ****f;
    double ****f_new;
    double ***rho;
    double ****u;

} lbm_params_3d;

void initialize_3d(lbm_params_3d lbm);
double equilibrium_3d(lbm_params_3d lbm, int k, double rho, double ux, double uy, double uz);
void collide_and_stream_3d(lbm_params_3d lbm);
void apply_boundary_conditions_3d_box(lbm_params_3d lbm);
void swap_distributions_3d(lbm_params_3d lbm);
void lbm_3d_step(lbm_params_3d lbm);
void lbm_free_3d(lbm_params_3d lbm);
