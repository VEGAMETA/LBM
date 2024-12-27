#pragma once
#include "utils.h"

typedef struct lbm_parameters_2d
{
    const int NX;
    const int NY;
    int D;
    const int Q;

    const double *w;

    const int *cx;
    const int *cy;

    const int *opp;
    const int *opp_x;
    const int *opp_y;

    const double EQ_A;
    const double EQ_B;
    const double EQ_C;
    const double EQ_D;

    double tau;

    double ***f;
    double ***f_new;
    double **rho;
    double ***u;
} lbm_params_2d;

typedef struct box{
    int x0;
    int x1;
    int y0;
    int y1;
} box;

void initialize_2d(lbm_params_2d *lbm);
double equilibrium_2d(lbm_params_2d lbm, int k, double rho, double ux, double uy);
void collide_and_stream_2d(lbm_params_2d lbm);
void boundary(lbm_params_2d *lbm, int x0, int x1, int y0, int y1);
void apply_boundary_conditions_2d(lbm_params_2d *lbm);
void swap_distributions_2d(lbm_params_2d lbm);
void lbm_2d_step(lbm_params_2d lbm);
void lbm_free_2d(lbm_params_2d lbm);