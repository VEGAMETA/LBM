#pragma once
#include <math.h>
#include <stdio.h>
#include <raylib.h>
#include "../lbm/lbm.h"

// Constants
#define GRID_SIZE_X 10
#define GRID_SIZE_Y 10
#define GRID_SIZE_Z 10

void generate_3d_data(float rho[GRID_SIZE_X][GRID_SIZE_Y][GRID_SIZE_Z]);
void generate_2d_data(float rho[GRID_SIZE_X][GRID_SIZE_Y]);

void visualize_3d(lbm_params_3d lbm);
void visualize_2d(lbm_params_2d lbm);