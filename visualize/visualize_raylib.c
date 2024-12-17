#include "visualize_raylib.h"

void visualize_2d(lbm_params_2d lbm) {
    const int screenWidth = 1200;
    const int screenHeight = 800;
    InitWindow(screenWidth, screenHeight, "2D Grid Visualization");

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 0.0f, -1000.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 95.0f;

    // Main loop
    while (!WindowShouldClose()) {
        lbm_2d_step(lbm);
        printf("rho[0][0]: %f\n", lbm.rho[0][0]); // ??
        // Draw
        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode3D(camera);
        // // Draw grid and data points
        #pragma omp parallel for collapse(2)
        for (int x = 0; x < lbm.NX; x++) {
            for (int y = 0; y < lbm.NY; y++) {
                Color color = ColorFromHSV((lbm.rho[x][y] - 1.) * 50.f+150.0f, 1.0f, 1.0f);
                DrawRectangle(4*(x-lbm.NX/2), 4*(y-lbm.NY/2), 5, 5, color);
            }
        }
        EndMode3D();
        DrawFPS(10, 10);
        EndDrawing();
    }

    CloseWindow();
};

void visualize_3d(lbm_params_3d lbm) {
    // Initialize the window
    const int screenWidth = 800;
    const int screenHeight = 600;
    InitWindow(screenWidth, screenHeight, "3D Grid Visualization");

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 20.0f, 20.0f, 20.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;

    // Data array
    static float rho[GRID_SIZE_X][GRID_SIZE_Y][GRID_SIZE_Z];
    //generate_3d_data(rho);

    // Main loop
    SetTargetFPS(60);
    while (!WindowShouldClose()) {




        // Update camera
        UpdateCamera(&camera, CAMERA_ORBITAL);

        // Draw
        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

        // // Draw grid and data points
        // for (int x = 0; x < GRID_SIZE_X; x++) {
        //     for (int y = 0; y < GRID_SIZE_Y; y++) {
        //         for (int z = 0; z < GRID_SIZE_Z; z++) {
        //             float value = rho[x][y][z];
        //             Color color = ColorFromHSV((value + 1.0f) * 180.0f, 1.0f, 1.0f);
        //             DrawCube((Vector3){ x, y, z }, 0.5f, 0.5f, 0.5f, color);
        //         }
        //     }
        // }

        EndMode3D();
        DrawFPS(10, 10);
        EndDrawing();
    }

    CloseWindow();
}
