#pragma once

double ****create_4d_array(int d1, int d2, int d3, int d4);
double ***create_3d_array(int d1, int d2, int d3);
double **create_2d_array(int d1, int d2);
double *create_1d_array(int d1);
void free_4d_array(double ****array, int d1, int d2, int d3);
void free_3d_array(double ***array, int d1, int d2);
void free_2d_array(double **array, int d1);
void free_1d_array(double *array);

typedef struct lbm_parameters
{
    const int NX;
    const int NY;
    const int NZ;
    int D;
    int Q;

    const double *w;

    const int *cx;
    const int *cy;
    const int *cz;

    const double EQ_A;
    const double EQ_B;
    const double EQ_C;
    const double EQ_D;

    double tau;
} lbm_params;

void initialize_3d(lbm_params lbm, double ****f, double ***rho, double ****u);
void initialize_2d(lbm_params lbm, double ***f, double **rho, double ***u);
void initialize_1d(lbm_params lbm, double **f, double *rho, double **u);

double equilibrium_3d(lbm_params lbm, int k, double rho, double ux, double uy, double uz);
double equilibrium_2d(lbm_params lbm, int k, double rho, double ux, double uy);
double equilibrium_1d(lbm_params lbm, int k, double rho, double ux);

void apply_boundary_conditions_3d_box(lbm_params lbm, double ****f);
void apply_boundary_conditions_2d_sandwich(lbm_params lbm, double ***f);
void apply_boundary_conditions_1d_ray(lbm_params lbm, double **f);

void collide_and_stream_3d(lbm_params lbm, double ****f, double ****f_new, double ***rho, double ****u);
void collide_and_stream_2d(lbm_params lbm, double ***f, double ***f_new, double **rho, double ***u);
void collide_and_stream_1d(lbm_params lbm, double **f, double **f_new, double *rho, double **u);

void swap_distributions_3d(lbm_params lbm, double ****f, double ****f_new);
void swap_distributions_2d(lbm_params lbm, double ***f, double ***f_new);
void swap_distributions_1d(lbm_params lbm, double **f, double **f_new);

void write_to_file_3d(lbm_params lbm, double ***rho);
void write_to_file_2d(lbm_params lbm, double **rho);
void write_to_file_1d(lbm_params lbm, double *rho);

void lbm_3d_loop(lbm_params lbm, int steps);
void lbm_2d_loop(lbm_params lbm, int steps);
void lbm_1d_loop(lbm_params lbm, int steps);
