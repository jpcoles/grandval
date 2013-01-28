#ifndef IO_H
#define IO_H

int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols);
void capture(double radius, struct particle *P, size_t NP, struct image *image, int clear, int rgb[3]);
void capture_massive(double radius, struct massive_particle *P, size_t NP, struct image *image, int clear, int rgb[3]);
void save_snapshot(int step, struct image *image);

#endif
