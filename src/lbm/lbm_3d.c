#include "lbm/lbm_3d.h"

void initialize_3d(lbm_params_3d lbm) {
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
}

double equilibrium_3d(lbm_params_3d lbm, int k, double rho, double ux, double uy, double uz) {
    double cu = lbm.cx[k] * ux + lbm.cy[k] * uy + lbm.cz[k] * uz;
    double usq = ux * ux + uy * uy + uz * uz;
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

void swap_distributions_3d(lbm_params_3d lbm){
    #pragma omp parallel for collapse(3)
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

void lbm_3d_step(lbm_params_3d lbm) {
    collide_and_stream_3d(lbm);
    apply_boundary_conditions_3d_box(lbm);
    swap_distributions_3d(lbm);
}

void lbm_free_3d(lbm_params_3d lbm){
    free_3d_array(lbm.rho, lbm.NX, lbm.NY);
    free_4d_array(lbm.f, lbm.NX, lbm.NY, lbm.NZ);
    free_4d_array(lbm.f_new, lbm.NX, lbm.NY, lbm.NZ);
    free_4d_array(lbm.u, lbm.NX, lbm.NY, lbm.NZ);   
}
