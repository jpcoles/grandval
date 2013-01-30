#ifndef IO_H
#define IO_H

#include "grandval.h"

struct io
{
    char *name;
    char *format;
    int overwrite;
};

void show_image_formats();
int check_image_format(char *f);

int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols);
void capture(dist_t radius, struct particle *P, size_t NP, struct image *image, int clear, int rgb[3]);
void capture_massive(dist_t radius, struct massive_particle *P, size_t NP, struct image *image, int clear, int rgb[3]);
void save_image(int step, struct image *image);

void show_snapshot_formats();
int check_snapshot_format(char *f);
int save_snapshot(int step, struct io *io, struct particle *P, size_t NP);

#endif

