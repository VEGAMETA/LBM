#include "lbm.h"

double ****create_4d_array(int d1, int d2, int d3, int d4) {
    double ****array = malloc(d1 * sizeof(double ***));
    for (int i = 0; i < d1; i++) {
        array[i] = malloc(d2 * sizeof(double **));
        for (int j = 0; j < d2; j++) {
            array[i][j] = malloc(d3 * sizeof(double *));
            for (int k = 0; k < d3; k++) {
                array[i][j][k] = malloc(d4 * sizeof(double));
            }
        }
    }
    return array;
}

double ***create_3d_array(int d1, int d2, int d3) {
    double ***array = malloc(d1 * sizeof(double **));
    for (int i = 0; i < d1; i++) {
        array[i] = malloc(d2 * sizeof(double *));
        for (int j = 0; j < d2; j++) {
            array[i][j] = malloc(d3 * sizeof(double));
        }
    }
    return array;
}

double **create_2d_array(int d1, int d2) {
    double **array = malloc(d1 * sizeof(double *));
    for (int i = 0; i < d1; i++) {
        array[i] = malloc(d2 * sizeof(double));
    }
    return array;
}

double *create_1d_array(int d1) {
    double *array = malloc(d1 * sizeof(double));
    return array;
}


void free_4d_array(double ****array, int d1, int d2, int d3) {
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            for (int k = 0; k < d3; k++) {
                free(array[i][j][k]);
            }
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

void free_3d_array(double ***array, int d1, int d2) {
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

void free_2d_array(double **array, int d1) {
    for (int i = 0; i < d1; i++) {
        free(array[i]);
    }
    free(array);
}

void free_1d_array(double *array) {
    free(array);
}


lbm_params_3d initialize_3d(lbm_params_3d lbm) {
    double ****f = create_4d_array(lbm.NX, lbm.NY, lbm.NZ, lbm.Q);
    double ****f_new = create_4d_array(lbm.NX, lbm.NY, lbm.NZ, lbm.Q);
    double ***rho = create_3d_array(lbm.NX, lbm.NY, lbm.NZ);
    double ****u = create_4d_array(lbm.NX, lbm.NY, lbm.NZ, lbm.D); 

    lbm.f = f;
    lbm.f_new = f_new;
    lbm.rho = rho;
    lbm.u = u;

    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                lbm.rho[x][y][z] = 1.0; // Initial density
                if (x*x + y*y > 200.0) lbm.rho[x][y][z] += 1.0;
                lbm.u[x][y][z][0] = 0.0;
                lbm.u[x][y][z][1] = 0.0;
                lbm.u[x][y][z][2] = 0.0;

                // Initialize equilibrium distribution
                for (int k = 0; k < lbm.Q; k++) {
                    lbm.f[x][y][z][k] = equilibrium_3d(lbm, k, lbm.rho[x][y][z], lbm.u[x][y][z][0], lbm.u[x][y][z][1], lbm.u[x][y][z][2]);
                }
            }
        }
    }
    return lbm;
}

void initialize_2d(lbm_params_2d *lbm) {
    double ***f = create_3d_array(lbm->NX, lbm->NY,lbm->Q);
    double ***f_new = create_3d_array(lbm->NX, lbm->NY, lbm->Q);
    // Macroscopic quantities
    double **rho = create_2d_array(lbm->NX, lbm->NY);
    double ***u = create_3d_array(lbm->NX, lbm->NY, lbm->D); 
    
    lbm->f = f;
    lbm->f_new = f_new;
    lbm->rho = rho;
    lbm->u = u;

    for (int x = 0; x < lbm->NX; x++) {
        for (int y = 0; y < lbm->NY; y++) {
            lbm->rho[x][y] = 1.0;
            lbm->rho[x][y] += 100*(((x-100)*(x-100) +  (y-250)*(y-250)) < 24000) || (x > 450 && y > 120 && y < 180); // Initial density
            lbm->rho[x][y] *= 1 + !(250 > x && x > 150 && 395 > y && y > 220);
            //if () lbm->rho[x][y] = 1.0;
            lbm->u[x][y][0] = 0.; // Initial x-velocity
            lbm->u[x][y][1] = 0.; // Initial y-velocity

            // Initialize equilibrium distribution
            for (int k = 0; k < lbm->Q; k++) {
                lbm->f[x][y][k] = equilibrium_2d(*lbm, k, lbm->rho[x][y], lbm->u[x][y][0], lbm->u[x][y][1]);
            }
        }
    }
}

lbm_params_1d initialize_1d(lbm_params_1d lbm) {
    double **f = create_2d_array(lbm.NX, lbm.Q);
    double **f_new = create_2d_array(lbm.NX, lbm.Q);
    // Macroscopic quantities
    double *rho = create_1d_array(lbm.NX);
    double **u = create_2d_array(lbm.NX, lbm.D); 

    lbm.f = f;
    lbm.f_new = f_new;
    lbm.rho = rho;
    lbm.u = u;

    for (int x = 0; x < lbm.NX; x++) {
        lbm.rho[x] = 1.0; // Initial density
        lbm.u[x][0] = 0.0; // Initial velocity

        // Initialize equilibrium distribution
        for (int k = 0; k < lbm.Q; k++) {
            lbm.f[x][k] = equilibrium_1d(lbm, k, lbm.rho[x], lbm.u[x][0]);
        }
    }
    return lbm;
}


double equilibrium_3d(lbm_params_3d lbm, int k, double rho, double ux, double uy, double uz) {
    double cu = lbm.cx[k] * ux + lbm.cy[k] * uy + lbm.cz[k] * uz;
    double usq = ux * ux + uy * uy + uz * uz;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

double equilibrium_2d(lbm_params_2d lbm, int k, double rho, double ux, double uy) {
    double cu = lbm.cx[k] * ux + lbm.cy[k] * uy;
    double usq = ux * ux + uy * uy;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

double equilibrium_1d(lbm_params_1d lbm, int k, double rho, double ux) {
    double cu = lbm.cx[k] * ux;
    double usq = ux * ux;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

void apply_face_3d(lbm_params_3d lbm, int fixed_coord, int fixed_value, int direction, int size1, int size2, const int *c1, const int *c2) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            for (int k = 0; k < lbm.Q; k++) {
                if (direction * c2[k] > 0) { 
                    int opp_k;
                    if (k == 0 && lbm.Q % 2 == 0) opp_k = 0;
                    else if ((k % 2 == 0 && lbm.Q % 2 == 0) || (k % 2 == 1 && lbm.Q % 2 == 1)) opp_k = k + 1;
                    else opp_k = k - 1;
                    if (fixed_coord == 0) {
                        lbm.f[fixed_value][i][j][k] = lbm.f[fixed_value][i][j][opp_k];
                    } else if (fixed_coord == 1) {
                        lbm.f[i][fixed_value][j][k] = lbm.f[i][fixed_value][j][opp_k];
                    } else {
                        lbm.f[i][j][fixed_value][k] = lbm.f[i][j][fixed_value][opp_k];
                    }
                }
            }
        }
    }
}

void apply_boundary_conditions_3d_box(lbm_params_3d lbm) {
    apply_face_3d(lbm, 2, 0, -1, lbm.NX, lbm.NY, lbm.cx, lbm.cz);            // Bottom face (Z=0)
    apply_face_3d(lbm, 2, lbm.NZ - 1, 1, lbm.NX, lbm.NY, lbm.cx, lbm.cz);    // Top face (Z=NZ-1)
    apply_face_3d(lbm, 0, 0, -1, lbm.NY, lbm.NZ, lbm.cx, lbm.cy);            // Left face (X=0)
    apply_face_3d(lbm, 0, lbm.NX - 1, 1, lbm.NY, lbm.NZ, lbm.cx, lbm.cy);    // Right face (X=NX-1)
    apply_face_3d(lbm, 1, 0, -1, lbm.NX, lbm.NZ, lbm.cx, lbm.cy);            // Front face (Y=0)
    apply_face_3d(lbm, 1, lbm.NY - 1, 1, lbm.NX, lbm.NZ, lbm.cx, lbm.cy);    // Back face (Y=NY-1)
}

void boundary(lbm_params_2d *lbm, int x, int y,const int *opposite, const int *c, int direction) {
    lbm->rho[x][y] = 0;
    lbm->u[x][y][0] = 0;
    lbm->u[x][y][1] = 0;
    double *temp = (double *)malloc(lbm->Q * sizeof(double));
    for (int q = 0; q < lbm->Q; q++) {
        temp[q] = lbm->f[x][y][q];
    }
    for (int q = 0; q < lbm->Q; q++) {
        lbm->f[x][y][q] = temp[opposite[q]];
        lbm->rho[x][y] += lbm->f[x][y][q];
        lbm->u[x][y][0] += lbm->cx[q] * lbm->f[x][y][q];
        lbm->u[x][y][1] += lbm->cy[q] * lbm->f[x][y][q];
    }
    free(temp);
    lbm->u[x][y][0] /= lbm->rho[x][y];
    lbm->u[x][y][1] /= lbm->rho[x][y];
}

void apply_inner_box_boundaries(lbm_params_2d *lbm, int x0, int x1, int y0, int y1) {
    for (int y = y0; y < y1; y++) {
            boundary(lbm, x0, y, lbm->opp_x, lbm->cx, -1); // right
            boundary(lbm, x1 - 1, y, lbm->opp_x, lbm->cx, 1); // left
    }
    for (int x = x0; x < x1; x++) {
            boundary(lbm, x, y0, lbm->opp_y, lbm->cy, -1); // bottom
            boundary(lbm, x, y1 - 1, lbm->opp_y, lbm->cy, 1); // top
    }
}

void apply_outer_box_boundaries(lbm_params_2d *lbm, int x0, int x1, int y0, int y1) {
    for (int y = y0; y < y1; y++) {
            boundary(lbm, x1 - 1, y, lbm->opp_x, lbm->cx, -1); // right
            boundary(lbm, x0, y, lbm->opp_x, lbm->cx, 1); // left
    }
    for (int x = x0; x < x1; x++) {
            boundary(lbm, x, y1 - 1, lbm->opp_y, lbm->cy, -1); // bottom
            boundary(lbm, x, y0, lbm->opp_y, lbm->cy, 1); // top
    }
}

void apply_boundary_conditions_2d(lbm_params_2d *lbm) {
    //apply_inner_box_boundaries(lbm, 150, 250, 220, 395);
    //apply_outer_box_boundaries(lbm, 150, 250, 220, 395);
    //apply_inner_box_boundaries(lbm, 0, lbm->NX, 0, lbm->NY);
    //apply_inner_box_boundaries(lbm, 150, 300, 220, 395);
}


void apply_boundary_conditions_1d_ray(lbm_params_1d lbm) {
    for (int k = 0; k < lbm.Q; k++) {
        lbm.f[0][k] = lbm.f[1][k];
    }
}


void collide_and_stream_3d(lbm_params_3d lbm){
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                // Compute macroscopic quantities
                lbm.rho[x][y][z] = 0.0;
                lbm.u[x][y][z][0] = 0.0;
                lbm.u[x][y][z][1] = 0.0;
                lbm.    u[x][y][z][2] = 0.0;
                for (int k = 0; k < lbm.Q; k++) {
                    lbm.rho[x][y][z] += lbm.f[x][y][z][k];
                    lbm.u[x][y][z][0] += lbm.f[x][y][z][k] * lbm.cx[k];
                    lbm.u[x][y][z][1] += lbm.f[x][y][z][k] * lbm.cy[k];
                    lbm.u[x][y][z][2] += lbm.f[x][y][z][k] * lbm.cz[k];
                }
                lbm.u[x][y][z][0] /= lbm.rho[x][y][z];
                lbm.u[x][y][z][1] /= lbm.rho[x][y][z];
                lbm.u[x][y][z][2] /= lbm.rho[x][y][z];

                //u[x][y][z][2] -= 9.8;

                // Collision step
                for (int k = 0; k < lbm.Q; k++) {
                    double feq = equilibrium_3d(lbm, k, lbm.rho[x][y][z], lbm.u[x][y][z][0], lbm.u[x][y][z][1], lbm.u[x][y][z][2]);
                    lbm.f_new[x][y][z][k] = lbm.f[x][y][z][k] - (lbm.f[x][y][z][k] - feq) / lbm.tau;
                }
            }
        }
    }
    // Streaming step
    #pragma omp parallel for collapse(4)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                for (int k = 0; k < lbm.Q; k++) {
                    int xp = (x + lbm.cx[k] + lbm.NX) % lbm.NX;
                    int yp = (y + lbm.cy[k] + lbm.NY) % lbm.NY;
                    int zp = (z + lbm.cz[k] + lbm.NZ) % lbm.NZ;
                    lbm.f_new[xp][yp][zp][k] = lbm.f[x][y][z][k];
                }
            }
        }
    }
}

void collide_and_stream_2d(lbm_params_2d lbm) {
    #pragma omp parallel for collapse(2)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            // Compute macroscopic quantities
            lbm.rho[x][y] = 0.0;
            lbm.u[x][y][0] = 0.0;
            lbm.u[x][y][1] = 0.0;
            for (int k = 0; k < lbm.Q; k++) {
                lbm.rho[x][y] += lbm.f[x][y][k];
                lbm.u[x][y][0] += lbm.f[x][y][k] * lbm.cx[k];
                lbm.u[x][y][1] += lbm.f[x][y][k] * lbm.cy[k];
            }
            if (lbm.rho[x][y] == 0.0) {
                lbm.u[x][y][0] /= 1e-10;
                lbm.u[x][y][1] /= 1e-10;
            }
            else {
                lbm.u[x][y][0] /= lbm.rho[x][y];
                lbm.u[x][y][1] /= lbm.rho[x][y];
            }

            // Collision step
            for (int k = 0; k < lbm.Q; k++) {
                double feq = equilibrium_2d(lbm, k, lbm.rho[x][y], lbm.u[x][y][0], lbm.u[x][y][1]);
                lbm.f[x][y][k] += (feq - lbm.f[x][y][k]) / lbm.tau;
            }
        }
    }

    // Streaming step
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int k = 0; k < lbm.Q; k++) {
                int xp = x + lbm.cx[k];
                int yp = y + lbm.cy[k];

                if (yp == -1 || yp == lbm.NY || xp == -1 || xp == lbm.NX) 
                    lbm.f_new[x][y][lbm.opp[k]] = lbm.f[x][y][k];
                else lbm.f_new[xp][yp][k] = lbm.f[x][y][k];
            }
        }
    }
}

void collide_and_stream_1d(lbm_params_1d lbm)
{
    for (int x = 0; x < lbm.NX; x++) {
        // Compute macroscopic quantities
        lbm.rho[x] = 0.0;
        lbm.u[x][0] = 0.0;
        for (int k = 0; k < 19; k++) {
            lbm.rho[x] += lbm.f[x][k];
            lbm.u[x][0] += lbm.f[x][k] * lbm.cx[k];
        }
        lbm.u[x][0] /= lbm.rho[x];

        // Collision step
        for (int k = 0; k < lbm.Q; k++) {
            double feq = equilibrium_1d(lbm, k, lbm.rho[x], lbm.u[x][0]);
            lbm.f_new[x][k] = lbm.f[x][k] - (lbm.f[x][k] - feq) / lbm.tau;
        }
    }

    // Streaming step
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            int xp = (x + lbm.cx[k] + lbm.NX) % lbm.NX;
            lbm.f_new[xp][k] = lbm.f[x][k];
        }       
    }
}


void swap_distributions_3d(lbm_params_3d lbm){
    #pragma omp parallel collapse(4)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                for (int k = 0; k < lbm.Q; k++) {
                    lbm.f[x][y][z][k] = lbm.f_new[x][y][z][k];
                }
            }
        }
    }
}

void swap_distributions_2d(lbm_params_2d lbm){
    #pragma omp parallel collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int k = 0; k < lbm.Q; k++) {
                lbm.f[x][y][k] = lbm.f_new[x][y][k];
            }
        }
    }
}

void swap_distributions_1d(lbm_params_1d lbm){
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            lbm.f[x][k] = lbm.f_new[x][k];
        }
    }
}


void write_to_file_3d(lbm_params_3d lbm) {
    FILE *file = fopen("density.dat", "w");
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                fprintf(file, "%d %d %d %f\n", x, y, z, lbm.rho[x][y][z]);
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_to_file_2d(lbm_params_2d lbm) {
    FILE *file = fopen("density.dat", "w");
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            fprintf(file, "%d %d %f\n", x, y, lbm.rho[x][y]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_to_file_1d(lbm_params_1d lbm) {
    FILE *file = fopen("density.dat", "w");
    for (int x = 0; x < lbm.NX; x++) {
            fprintf(file, "%d %f\n", x, lbm.rho[x]);
    }
    fclose(file);
}


void lbm_3d_step(lbm_params_3d lbm) {
    collide_and_stream_3d(lbm);
    apply_boundary_conditions_3d_box(lbm);
    swap_distributions_3d(lbm);
}

void lbm_2d_step(lbm_params_2d lbm) {
    collide_and_stream_2d(lbm);
    apply_boundary_conditions_2d(&lbm);
    swap_distributions_2d(lbm);
}


void lbm_1d_step(lbm_params_1d lbm) {
    collide_and_stream_1d(lbm);
    apply_boundary_conditions_1d_ray(lbm);
    swap_distributions_1d(lbm);
}

void lbm_free_3d(lbm_params_3d lbm){
    free_3d_array(lbm.rho, lbm.NX, lbm.NY);
    free_4d_array(lbm.f, lbm.NX, lbm.NY, lbm.NZ);
    free_4d_array(lbm.f_new, lbm.NX, lbm.NY, lbm.NZ);
    free_4d_array(lbm.u, lbm.NX, lbm.NY, lbm.NZ);   
}

void lbm_free_2d(lbm_params_2d lbm) {
    free_2d_array(lbm.rho, lbm.NX);
    free_3d_array(lbm.f, lbm.NX, lbm.NY);
    free_3d_array(lbm.f_new, lbm.NX, lbm.NY);
    free_3d_array(lbm.u, lbm.NX, lbm.NY);
}

void lbm_free_1d(lbm_params_1d lbm){
    free_1d_array(lbm.rho);
    free_2d_array(lbm.f, lbm.NX);
    free_2d_array(lbm.f_new, lbm.NX);
    free_2d_array(lbm.u, lbm.NX);
}