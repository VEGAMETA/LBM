#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#define _USE_MATH_DEFINES

double ****create_4d_array(int d1, int d2, int d3, int d4);
double ***create_3d_array(int d1, int d2, int d3);
double **create_2d_array(int d1, int d2);
double *create_1d_array(int d1);
void free_4d_array(double ****array, int d1, int d2, int d3);
void free_3d_array(double ***array, int d1, int d2);
void free_2d_array(double **array, int d1);
void free_1d_array(double *array);

double clamp(double x, double min, double max);