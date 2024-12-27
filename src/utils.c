#include "utils.h"

double ****create_4d_array(int d1, int d2, int d3, int d4) {
    double ****array = (double ****)malloc(d1 * sizeof(double ***));
    for (int i = 0; i < d1; i++) array[i] = create_3d_array(d2, d3, d4);
    return array;
}

double ***create_3d_array(int d1, int d2, int d3) {
    double ***array = (double ***)malloc(d1 * sizeof(double **));
    for (int i = 0; i < d1; i++) array[i] = create_2d_array(d2, d3);
    return array;
}

double **create_2d_array(int d1, int d2) {
    double **array = (double **)malloc(d1 * sizeof(double *));
    for (int i = 0; i < d1; i++) array[i] = create_1d_array(d2);
    return array;
}

double *create_1d_array(int d1) {
    return (double *)malloc(d1 * sizeof(double));
}


void free_4d_array(double ****array, int d1, int d2, int d3) {
    for (int i = 0; i < d1; i++) free_3d_array(array[i], d2, d3);
    free(array);
}

void free_3d_array(double ***array, int d1, int d2) {
    for (int i = 0; i < d1; i++) free_2d_array(array[i], d2);
    free(array);
}

void free_2d_array(double **array, int d1) {
    for (int i = 0; i < d1; i++) free_1d_array(array[i]);
    free(array);
}

void free_1d_array(double *array) {
    free(array);
}

double clamp(double value, double min, double max) {
    return fmin(fmax(value, min), max);
}