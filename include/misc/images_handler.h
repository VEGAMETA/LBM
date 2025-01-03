#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jpeglib.h>
#include <png.h>
#include <errno.h>
#include "lbm/utils.h"

void process_jpeg(const char *filename, int ***lattice_array, int *height, int *width);
void process_png(const char *filename, int ***lattice_array, int *height, int *width);
void process_image(const char *filename, int ***lattice_array, int *height, int *width);