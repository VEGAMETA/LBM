#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _USE_MATH_DEFINES

double ****create_4d_array(int d1, int d2, int d3, int d4);
double ***create_3d_array(int d1, int d2, int d3);
double **create_2d_array(int d1, int d2);
double *create_1d_array(int d1);
void free_4d_array(double ****array, int d1, int d2, int d3);
void free_3d_array(double ***array, int d1, int d2);
void free_2d_array(double **array, int d1);
void free_1d_array(double *array);

typedef enum WALL_TYPE
{
    RIGHT = 1,
    LEFT = -1,
    TOP = -1,
    BOTTOM = 1,
    FRONT = -1,
    BACK = 1
} wall_type;


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

typedef struct lbm_parameters_1d
{
    const int NX;
    int D;
    const int Q;

    const double *w;

    const int *cx;
    const double EQ_A;
    const double EQ_B;
    const double EQ_C;
    const double EQ_D;

    double tau;

    double **f;
    double **f_new;
    double *rho;
    double **u;
} lbm_params_1d;

lbm_params_3d initialize_3d(lbm_params_3d lbm);
void initialize_2d(lbm_params_2d *lbm);
lbm_params_1d initialize_1d(lbm_params_1d lbm);

double equilibrium_3d(lbm_params_3d lbm, int k, double rho, double ux, double uy, double uz);
double equilibrium_2d(lbm_params_2d lbm, int k, double rho, double ux, double uy);
double equilibrium_1d(lbm_params_1d lbm, int k, double rho, double ux);

void boundary(lbm_params_2d *lbm, int x, int y, const int *opposite, const int *c, int direction);
void apply_box_boundaries(lbm_params_2d *lbm, int x0, int x1, int y0, int y1);

void apply_boundary_conditions_3d_box(lbm_params_3d lbm);
void apply_boundary_conditions_2d(lbm_params_2d *lbm);
void apply_boundary_conditions_1d_ray(lbm_params_1d lbm);

void collide_and_stream_3d(lbm_params_3d lbm);
void collide_and_stream_2d(lbm_params_2d lbm);
void collide_and_stream_1d(lbm_params_1d lbm);

void swap_distributions_3d(lbm_params_3d lbm);
void swap_distributions_2d(lbm_params_2d lbm);
void swap_distributions_1d(lbm_params_1d lbm);

void write_to_file_3d(lbm_params_3d lbm);
void write_to_file_2d(lbm_params_2d lbm);
void write_to_file_1d(lbm_params_1d lbm);

void lbm_3d_step(lbm_params_3d lbm);
void lbm_2d_step(lbm_params_2d lbm);
void lbm_1d_step(lbm_params_1d lbm);

void lbm_free_2d(lbm_params_2d lbm);
void lbm_free_3d(lbm_params_3d lbm);
void lbm_free_1d(lbm_params_1d lbm);
