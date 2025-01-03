#include "visualization/visualize.h"

void visualize_1d(lbm_params_1d lbm) {
    const int screenWidth = 1200;
    const int screenHeight = 800;
    InitWindow(screenWidth, screenHeight, "1D LBM Visualization");
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 0.0f, -1000.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 95.0f;
    SetTargetFPS(5);
    while (!WindowShouldClose()) {
        lbm_1d_step(lbm);
        printf_s("lbm.rho: %f\n", lbm.rho[lbm.NX - 2]);
        BeginDrawing(); 
        ClearBackground(BLACK);
        BeginMode3D(camera);
        #pragma omp parallel for
        for (int x = 0; x < lbm.NX; x++) {
            Color color = ColorFromHSV((lbm.rho[x] - 1.) * 100.f+100.0f, 1.0f, clamp(lbm.rho[x] - 1., 0.0f, 1.0f));
            DrawRectangle((-4*x) + lbm.NX * 2, 0, 4, 350, color);
        }
        EndMode3D();
        DrawFPS(10, 10);
        EndDrawing();
    }
}

void visualize_2d(lbm_params_2d lbm) {
    const int screenWidth = 1200;
    const int screenHeight = 800;
    InitWindow(screenWidth, screenHeight, "2D LBM Visualization");
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 0.0f, -1000.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 90.0f;
    //SetTargetFPS(80);
    while (!WindowShouldClose()) {
        lbm_2d_step(lbm);
        BeginDrawing();
        ClearBackground(BLACK);
        BeginMode3D(camera);
        for (int x = 0; x < lbm.NX; x++) {
            for (int y = 0; y < lbm.NY; y++) {
                Color color = ColorFromHSV((2 - lbm.rho[x][y]) * 180.f, 1.0f, clamp(lbm.rho[x][y], 0.0f, 1.0f));
                DrawRectangle((-4*x)+lbm.NX*2, (4*y)-lbm.NY*2, 5, 5, color);
            }
        }
        EndMode3D();
        DrawFPS(10, 10);
        EndDrawing();
    }
    CloseWindow();
};

void visualize_3d(lbm_params_3d lbm) {
    const int screenWidth = 800;
    const int screenHeight = 600;
    InitWindow(screenWidth, screenHeight, "3D Grid Visualization");
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 20.0f, 20.0f, 20.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        UpdateCamera(&camera, CAMERA_ORBITAL);
        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode3D(camera);

        EndMode3D();
        DrawFPS(10, 10);
        EndDrawing();
    }
    CloseWindow();
}
