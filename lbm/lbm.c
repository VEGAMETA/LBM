#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lbm.h"

#define _USE_MATH_DEFINES

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


void initialize_3d(lbm_params lbm, double ****f, double ***rho, double ****u) {
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                rho[x][y][z] = 1.0; // Initial density
                u[x][y][z][0] = 0.0;
                u[x][y][z][1] = 0.0;
                u[x][y][z][2] = 0.0;

                // Initialize equilibrium distribution
                for (int k = 0; k < lbm.Q; k++) {
                    f[x][y][z][k] = equilibrium_3d(lbm, k, rho[x][y][z], u[x][y][z][0], u[x][y][z][1], u[x][y][z][2]);
                }
            }
        }
    }
}

void initialize_2d(lbm_params lbm, double ***f, double **rho, double ***u) {
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            rho[x][y] = 1.0; // Initial density
            u[x][y][0] = 0.0; // Initial x-velocity
            u[x][y][1] = 0.0; // Initial y-velocity

            // Initialize equilibrium distribution
            for (int k = 0; k < lbm.Q; k++) {
                f[x][y][k] = equilibrium_2d(lbm, k, rho[x][y], u[x][y][0], u[x][y][1]);
            }
        }
    }
}

void initialize_1d(lbm_params lbm, double **f, double *rho, double **u) {
    for (int x = 0; x < lbm.NX; x++) {
        rho[x] = 1.0; // Initial density
        u[x][0] = 0.0; // Initial velocity

        // Initialize equilibrium distribution
        for (int k = 0; k < lbm.Q; k++) {
            f[x][k] = equilibrium_1d(lbm, k, rho[x], u[x][0]);
        }
    }
}


double equilibrium_3d(lbm_params lbm, int k, double rho, double ux, double uy, double uz) {
    double cu = lbm.cx[k] * ux + lbm.cy[k] * uy + lbm.cz[k] * uz;
    double usq = ux * ux + uy * uy + uz * uz;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

double equilibrium_2d(lbm_params lbm, int k, double rho, double ux, double uy) {
    double cu = lbm.cx[k] * ux + lbm.cy[k] * uy;
    double usq = ux * ux + uy * uy;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};

double equilibrium_1d(lbm_params lbm, int k, double rho, double ux) {
    double cu = lbm.cx[k] * ux;
    double usq = ux * ux;
    return lbm.w[k] * rho * (lbm.EQ_A + lbm.EQ_B * cu + lbm.EQ_C * cu * cu + lbm.EQ_D * usq);
};


void apply_boundary_conditions_3d_box(lbm_params lbm, double ****f) {
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            // Bottom face (Z=0)
            for (int k = 0; k < 19; k++) {
                if (lbm.cz[k] < 0) { // Opposite direction
                    int opp_k = 18 - k;
                    f[x][y][0][k] = f[x][y][0][opp_k];
                }
            }
            // Top face (Z=NZ-1)
            for (int k = 0; k < 19; k++) {
                if (lbm.cz[k] > 0) { // Opposite direction
                    int opp_k = 18 - k;
                    f[x][y][lbm.NZ - 1][k] = f[x][y][lbm.NZ - 1][opp_k];
                }
            }
        }
    }

    // Left and Right (X=0 and X=NX-1)
    for (int y = 0; y < lbm.NY; y++) {
        for (int z = 0; z < lbm.NZ; z++) {
            // Left face (X=0)
            for (int k = 0; k < 19; k++) {
                if (lbm.cx[k] < 0) { // Opposite direction
                    int opp_k = 18 - k;
                    f[0][y][z][k] = f[0][y][z][opp_k];
                }
            }
            // Right face (X=NX-1)
            for (int k = 0; k < 19; k++) {
                if (lbm.cx[k] > 0) { // Opposite direction
                    int opp_k = 18 - k;
                    f[lbm.NX - 1][y][z][k] = f[lbm.NX - 1][y][z][opp_k];
                }
            }
        }
    }

    // Front and Back (Y=0 and Y=NY-1)
    for (int x = 0; x < lbm.NX; x++) {
        for (int z = 0; z < lbm.NZ; z++) {
            // Front face (Y=0)
            for (int k = 0; k < 19; k++) {
                if (lbm.cy[k] < 0) { // Opposite direction
                    int opp_k = 18 - k;
                    f[x][0][z][k] = f[x][0][z][opp_k];
                }
            }
            // Back face (Y=NY-1)
            for (int k = 0; k < 19; k++) {
                if (lbm.cy[k] > 0) { // Opposite direction
                    int opp_k = 18 - k;
                    f[x][lbm.NY - 1][z][k] = f[x][lbm.NY - 1][z][opp_k];
                }
            }
        }
    }
};

void apply_boundary_conditions_2d_sandwich(lbm_params lbm, double ***f) {
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            // Bottom wall
            f[x][0][k] = f[x][1][k];
            // Top wall
            f[x][lbm.NY-1][k] = f[x][lbm.NY-2][k];
        }
    }
};

void apply_boundary_conditions_1d_ray(lbm_params lbm, double **f) {
    for (int k = 0; k < lbm.Q; k++) {
        f[0][k] = f[1][k];
    }
};


void collide_and_stream_3d(lbm_params lbm, double ****f, double ****f_new, double ***rho, double ****u){
    #pragma omp parallel for collapse(3)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                // Compute macroscopic quantities
                rho[x][y][z] = 0.0;
                u[x][y][z][0] = 0.0;
                u[x][y][z][1] = 0.0;
                u[x][y][z][2] = 0.0;
                for (int k = 0; k < lbm.Q; k++) {
                    rho[x][y][z] += f[x][y][z][k];
                    u[x][y][z][0] += f[x][y][z][k] * lbm.cx[k];
                    u[x][y][z][1] += f[x][y][z][k] * lbm.cy[k];
                    u[x][y][z][2] += f[x][y][z][k] * lbm.cz[k];
                }
                u[x][y][z][0] /= rho[x][y][z];
                u[x][y][z][1] /= rho[x][y][z];
                u[x][y][z][2] /= rho[x][y][z];

                // Collision step
                for (int k = 0; k < lbm.Q; k++) {
                    double feq = equilibrium_3d(lbm, k, rho[x][y][z], u[x][y][z][0], u[x][y][z][1], u[x][y][z][2]);
                    f_new[x][y][z][k] = f[x][y][z][k] - (f[x][y][z][k] - feq) / lbm.tau;
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
                    f_new[xp][yp][zp][k] = f[x][y][z][k];
                }
            }
        }
    }
}

void collide_and_stream_2d(lbm_params lbm, double ***f, double ***f_new, double **rho, double ***u) {
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            // Compute macroscopic quantities
            rho[x][y] = 0.0;
            u[x][y][0] = 0.0;
            u[x][y][1] = 0.0;
            for (int k = 0; k < lbm.Q; k++) {
                rho[x][y] += f[x][y][k];
                u[x][y][0] += f[x][y][k] * lbm.cx[k];
                u[x][y][1] += f[x][y][k] * lbm.cy[k];
            }
            u[x][y][0] /= rho[x][y];
            u[x][y][1] /= rho[x][y];

            // Collision step
            for (int k = 0; k < lbm.Q; k++) {
                double feq = equilibrium_2d(lbm, k, rho[x][y], u[x][y][0], u[x][y][1]);
                f_new[x][y][k] = f[x][y][k] - (f[x][y][k] - feq) / lbm.tau;
            }
        }
    }

    // Streaming step
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int k = 0; k < lbm.Q; k++) {
                int xp = (x + lbm.cx[k] + lbm.NX) % lbm.NX;
                int yp = (y + lbm.cy[k] + lbm.NY) % lbm.NY;
                f_new[xp][yp][k] = f[x][y][k];
            }
        }
    }
}

void collide_and_stream_1d(lbm_params lbm, double **f, double **f_new, double *rho, double **u)
{
    for (int x = 0; x < lbm.NX; x++) {
        // Compute macroscopic quantities
        rho[x] = 0.0;
        u[x][0] = 0.0;
        for (int k = 0; k < 19; k++) {
            rho[x] += f[x][k];
            u[x][0] += f[x][k] * lbm.cx[k];
        }
        u[x][0] /= rho[x];

        // Collision step
        for (int k = 0; k < lbm.Q; k++) {
            double feq = equilibrium_1d(lbm, k, rho[x], u[x][0]);
            f_new[x][k] = f[x][k] - (f[x][k] - feq) / lbm.tau;
        }
    }

    // Streaming step
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            int xp = (x + lbm.cx[k] + lbm.NX) % lbm.NX;
            f_new[xp][k] = f[x][k];
        }       
    }
}


void swap_distributions_3d(lbm_params lbm, double ****f, double ****f_new){
    #pragma omp parallel collapse(4)
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                for (int k = 0; k < lbm.Q; k++) {
                    f[x][y][z][k] = f_new[x][y][z][k];
                }
            }
        }
    }
}

void swap_distributions_2d(lbm_params lbm, double ***f, double ***f_new){
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int k = 0; k < lbm.Q; k++) {
                f[x][y][k] = f_new[x][y][k];
            }
        }
    }
}

void swap_distributions_1d(lbm_params lbm, double **f, double **f_new){
    for (int x = 0; x < lbm.NX; x++) {
        for (int k = 0; k < lbm.Q; k++) {
            f[x][k] = f_new[x][k];
        }
    }
}


void write_to_file_3d(lbm_params lbm, double ***rho) {
    FILE *file = fopen("density.dat", "w");
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            for (int z = 0; z < lbm.NZ; z++) {
                fprintf(file, "%d %d %d %f\n", x, y, z, rho[x][y][z]);
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_to_file_2d(lbm_params lbm, double **rho) {
    FILE *file = fopen("density.dat", "w");
    for (int x = 0; x < lbm.NX; x++) {
        for (int y = 0; y < lbm.NY; y++) {
            fprintf(file, "%d %d %f\n", x, y, rho[x][y]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_to_file_1d(lbm_params lbm, double *rho) {
    FILE *file = fopen("density.dat", "w");
    for (int x = 0; x < lbm.NX; x++) {
            fprintf(file, "%d %f\n", x, rho[x]);
    }
    fclose(file);
}


void lbm_3d_loop(lbm_params lbm, int steps) {
    // Distribution functions
    double ****f = create_4d_array(lbm.NX, lbm.NY, lbm.NZ, lbm.Q);
    double ****f_new = create_4d_array(lbm.NX, lbm.NY, lbm.NZ, lbm.Q);
    // Macroscopic quantities
    double ***rho = create_3d_array(lbm.NX, lbm.NY, lbm.NZ);
    double ****u = create_4d_array(lbm.NX, lbm.NY, lbm.NZ, lbm.D); 

    initialize_3d(lbm, f, rho, u);

    for (int step = 0; step < steps; step++) {
        collide_and_stream_3d(lbm, f, f_new, rho, u);
        apply_boundary_conditions_3d_box(lbm, f_new);
        swap_distributions_3d(lbm, f, f_new);
    }

    write_to_file_3d(lbm, rho);

    free_3d_array(rho, lbm.NX, lbm.NY);
    free_4d_array(f, lbm.NX, lbm.NY, lbm.NZ);
    free_4d_array(f_new, lbm.NX, lbm.NY, lbm.NZ);
    free_4d_array(u, lbm.NX, lbm.NY, lbm.NZ);   
}

void lbm_2d_loop(lbm_params lbm, int steps) {
    double ***f = create_3d_array(lbm.NX, lbm.NY,lbm.Q);
    double ***f_new = create_3d_array(lbm.NX, lbm.NY, lbm.Q);
    // Macroscopic quantities
    double **rho = create_2d_array(lbm.NX, lbm.NY);
    double ***u = create_3d_array(lbm.NX, lbm.NY, lbm.D); 
    initialize_2d(lbm, f, rho, u);
   // Main time-stepping loop
    for (int step = 0; step < steps; step++) {
        collide_and_stream_2d(lbm, f, f_new, rho, u);
        apply_boundary_conditions_2d_sandwich(lbm, f_new);
        swap_distributions_2d(lbm, f, f_new);
    }
    // Output final density for visualization
    write_to_file_2d(lbm, rho);

    free_2d_array(rho, lbm.NX);
    free_3d_array(f, lbm.NX, lbm.NY);
    free_3d_array(f_new, lbm.NX, lbm.NY);
    free_3d_array(u, lbm.NX, lbm.NY);
}

void lbm_1d_loop(lbm_params lbm, int steps) {
    double **f = create_2d_array(lbm.NX, lbm.Q);
    double **f_new = create_2d_array(lbm.NX, lbm.Q);
    // Macroscopic quantities
    double *rho = create_1d_array(lbm.NX);
    double **u = create_2d_array(lbm.NX, lbm.D); 
    initialize_1d(lbm, f, rho, u);
   // Main time-stepping loop
    for (int step = 0; step < steps; step++) {
        collide_and_stream_1d(lbm, f, f_new, rho, u);
        apply_boundary_conditions_1d_ray(lbm, f_new);
        swap_distributions_1d(lbm, f, f_new);
    }

    // Output final density for visualization
    write_to_file_1d(lbm, rho);

    free_1d_array(rho);
    free_2d_array(f, lbm.NX);
    free_2d_array(f_new, lbm.NX);
    free_2d_array(u, lbm.NX);
}