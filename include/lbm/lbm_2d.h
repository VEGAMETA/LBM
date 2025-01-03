#pragma once
#include "lbm/utils.h"

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

    double tau;
    
    int **map;

    double ***f;
    double ***f_new;
    double **rho;
    double ***u;
} lbm_params_2d;

typedef struct{
    int x0;
    int x1;
    int y0;
    int y1;
} box;

typedef struct{
    int x;
    int y;
    double r;
} circle;

void initialize_2d(lbm_params_2d *lbm);
void map_initialize_2d(lbm_params_2d lbm);
void base_initialize_2d(lbm_params_2d lbm);
void lbm_2d_step(lbm_params_2d lbm);
void lbm_free_2d(lbm_params_2d lbm);
double equilibrium_2d(lbm_params_2d lbm, int q, double rho, double ux, double uy);
