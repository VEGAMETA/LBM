#include "lbm/utils.h"

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

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>

double get_time_in_seconds() {
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / frequency.QuadPart;
}
#else
#include <time.h>

double get_time_in_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}
#endif

double time_sin() {
    return sin(2. * M_PI * get_time_in_seconds());
}

double time_cos() {
    return cos(2. * M_PI * get_time_in_seconds());
}