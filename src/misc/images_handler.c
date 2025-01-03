#include "misc/images_handler.h"

void process_jpeg(const char *filename, int ***lattice_array, int *height, int *width) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return perror("File opening failed");
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, fp);
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);
    *width = cinfo.output_width;
    *height = cinfo.output_height;
    int num_channels = cinfo.output_components;
    unsigned long row_stride = *width * num_channels;
    unsigned char *buffer = (unsigned char *)malloc(row_stride);

    *lattice_array = (int **)malloc(sizeof(int *) * (*width));
    for (int i = 0; i < *width; i++) (*lattice_array)[i] = (int *)malloc(sizeof(int) * (*height));

    while (cinfo.output_scanline < cinfo.output_height) {
        jpeg_read_scanlines(&cinfo, &buffer, 1);
        for (unsigned int x = 0; x < *width; x++) {
            unsigned char *px = &buffer[x * num_channels];
            int y = *height - cinfo.output_scanline - 1;
            if (px[0] <= 100 && px[1] <= 100 && px[2] <= 100) (*lattice_array)[x][y] = ADIABATIC;
            else if (px[0] > 100 && px[1] <= 100 && px[2] <= 100) (*lattice_array)[x][y] = CONVECTIVE_OUT;
            else if (px[0] <= 100 && px[1] <= 100 && px[2] > 100) (*lattice_array)[x][y] = CONVECTIVE_IN;
            else (*lattice_array)[x][y] = SPACE;
        }
    }
    free(buffer);
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(fp);
}


void error_handler(const char *msg, FILE *fp, png_structp png, png_infop info) {
    if (errno != 0) perror(msg);
    else fprintf(stderr, "%s\n", msg);
    if (fp) fclose(fp);
    if (png && info) png_destroy_read_struct(&png, &info, NULL);
    else if (png) png_destroy_read_struct(&png, NULL, NULL);
}

void process_png(const char *filename, int ***lattice_array, int *height, int *width) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return error_handler("File opening failed", NULL, NULL, NULL);
    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) return error_handler("png_create_read_struct failed\n", fp, NULL, NULL);
    png_infop info = png_create_info_struct(png);
    if (!info) return error_handler("png_create_info_struct failed\n", fp, png, NULL);
    if (setjmp(png_jmpbuf(png))) return error_handler("Error during init_io\n", fp, png, info);
    png_init_io(png, fp);
    png_read_info(png, info);
    *width = png_get_image_width(png, info);
    *height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth = png_get_bit_depth(png, info);
    if (bit_depth == 16) png_set_strip_16(png);
    if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
    if (png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);
    png_read_update_info(png, info);
    png_bytep *row_pointers = (png_bytep *)malloc(sizeof(png_bytep) * (*height));
    for (int y = 0; y < *height; y++) row_pointers[y] = (png_byte *)malloc(png_get_rowbytes(png, info));
    png_read_image(png, row_pointers);

    *lattice_array = (int **)malloc(sizeof(int *) * (*width));
    for (int i = 0; i < *width; i++) (*lattice_array)[i] = (int *)malloc(sizeof(int) * (*height));

    for (int y = 0; y < *height; y++) {
        png_bytep row = row_pointers[y];
        for (int x = 0; x < *width; x++) {
            png_bytep px = &(row[x * 4]);
            if (px[0] <= 100 && px[1] <= 100 && px[2] <= 100) (*lattice_array)[x][*height - y] = ADIABATIC;
            else if (px[0] > 100 && px[1] <= 100 && px[2] <= 100) (*lattice_array)[x][*height - y] = CONVECTIVE_OUT;
            else if (px[0] <= 100 && px[1] <= 100 && px[2] > 100) (*lattice_array)[x][*height - y] = CONVECTIVE_IN;
            else (*lattice_array)[x][*height - y] = SPACE;
        }
    }
    for (int y = 0; y < *height; y++) free(row_pointers[y]);
    free(row_pointers);
    fclose(fp);
    png_destroy_read_struct(&png, &info, NULL);
}

void process_image(const char *filename, int ***lattice_array, int *height, int *width){
    const char *dot = strrchr(filename, '.');
    if (!dot || dot == filename) {
        fprintf(stderr, "Invalid file name\n");
        return;
    }
    dot = dot + 1;
    if (strcmp(dot, "jpg") == 0 || strcmp(dot, "jpeg") == 0) process_jpeg(filename, lattice_array, height, width);
    else if (strcmp(dot, "png") == 0) process_png(filename, lattice_array, height, width);
    else fprintf(stderr, "Not a supported image format\n");
}