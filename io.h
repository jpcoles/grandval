#ifndef IO_H
#define IO_H

#include "grandval.h"

struct io
{
    char *name;
    char *format;
    int overwrite;
};

#define MIN_HEADER_SIZE (1024)

struct binary_header
{
    char magic[8];
    uint16_t    major_version;
    uint16_t    minor_version;
    uint32_t    sizeof_header;
    uint64_t    creation_time;
    uint8_t     sizeof_pos_t;
    uint8_t     sizeof_vel_t;
    uint8_t     sizeof_energy_t;
    uint8_t     sizeof_id;
    uint64_t    simulation_step;
    double      simulation_time;
    uint64_t    output_index;
    uint64_t    Nparticles;
    uint64_t    sizeof_options;
    unsigned char options[0];
};

void show_image_formats();
int check_image_format(char *f);

int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols);
void capture(dist_t radius, struct particle *P, size_t NP, struct image *image, int clear, int rgb[3]);
void capture_massive(dist_t radius, struct massive_particle *P, size_t NP, struct image *image, int clear, int rgb[3]);
void save_image(int step, struct image *image);

void show_snapshot_formats();
int check_snapshot_format(char *f);
int save_snapshot(int step, struct io *io, struct binary_header *hdr, struct particle *P, size_t NP);

void create_io_header(struct binary_header **hdr, double sim_time, int sim_save, int sim_step, struct program_options *opts, struct program_options *default_opts, size_t NP);
void show_binary_snapshot(char *fname);


#endif

