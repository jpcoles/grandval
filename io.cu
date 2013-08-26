#include <endian.h>
#include <sys/time.h>
#include <errno.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <arpa/inet.h>
#include "grandval.h"

#include "io.h"
#include "util.h"
#include "color_ramp.h"
#include "options.h"

static struct
{
    char *name;
    char *desc;
} supported_image_formats[] = { 
    {"png", "Portable Network Graphic"}, 
    {NULL, NULL} 
};

static struct
{
    char *name;
    char *desc;
} supported_snapshot_formats[] = { 
    {"ascii", "ASCII Text. Space separated values. Header in comments."}, 
    {"csv", "ASCII Text. Comma separated values. Header in comments."}, 
    {"binary", "Binary format. Header is at least 1024 bytes long followed particle structures"}, 
    {NULL, NULL} 
};

void create_io_header(struct binary_header **hdr, double sim_time, int sim_save, int sim_step, struct program_options *opts, struct program_options *default_opts, size_t NP)
{
    struct binary_header *h = *hdr;

    if (h == NULL)
    {
        //char *opts_text = make_binary_options_text(opts, default_opts);
        //size_t sizeof_options = strlen(opts_text) + 1;

        //*hdr = (struct binary_header *)malloc(sizeof(**hdr) + sizeof_options*sizeof(*opts_text));
        *hdr = (struct binary_header *)malloc(sizeof(**hdr));
        h = *hdr;
        strncpy(h->magic, "GRANDVAL", 8);
        h->major_version = htole16(MAJOR_VERSION);
        h->minor_version = htole16(MINOR_VERSION);
        h->sizeof_header = htole32(MIN_HEADER_SIZE < sizeof(*h) ? sizeof(*h) : MIN_HEADER_SIZE);
        h->sizeof_pos_t  = (uint8_t)sizeof(pos_t);
        h->sizeof_vel_t  = (uint8_t)sizeof(vel_t);
        h->sizeof_energy_t=(uint8_t)sizeof(energy_t);
        h->sizeof_id     =  sizeof(int); 
        h->Nparticles = htole64(NP);

        h->sizeof_options = 0; //htonl(sizeof_options);
        //memcpy(h->options, opts_text, h->sizeof_options*sizeof(*opts_text));
    }

    h->creation_time   = htole64((uint64_t)time(NULL));
    h->simulation_step = htole64(sim_step);
    h->simulation_time = sim_time;
    h->output_index    = htole64(sim_save);
}

//==============================================================================
//                                  capture
//==============================================================================
void capture(dist_t radius, struct particle *P, size_t NP, struct image *image, int clear, int rgb[3])
{
    int32_t r, c;
    double scale = 100;

    if (clear) 
    {
        memset(image->image, 0, 3*image->nc*image->nr*sizeof(*(image->image)));
        memset(image->hist,  0, 1*image->nc*image->nr*sizeof(*(image->hist)));
    }

    //if (with_gradient)
        //scale = pow(((double)env->step) / env->opt.Nsteps, 3);

    size_t i;
    for (i=0; i < NP; i++)
    {
        c = ( P[i].x[0] + radius) / (2.0*radius) * image->nc;
        r = (-P[i].x[1] + radius) / (2.0*radius) * image->nr;
        if (!(0 <= r && r < image->nc)) continue;
        if (!(0 <= c && c < image->nr)) continue;
        image->hist[r*image->nc + c] += 1;
        //if (image->hist[r*image->nc + c] > scale)
            //scale = image->hist[r*image->nc + c];
    }

    i = 0;
    for (r=0; r < image->nr; r++)
    {
        for (c=0; c < image->nc; c++)
        {
            float r,g,b;

            r = g = b = 0;
            if (image->hist[i] > 0)
            {
                r = g = b = scale > 0 ? image->hist[i] / scale : 0;
                color_ramp_astro(&r,&g,&b);
            }

            image->image[3*i + 0] = r * 255;
            image->image[3*i + 1] = g * 255;
            image->image[3*i + 2] = b * 255;

#if 0
            image->image[3*i + 0] = scale > 0 ? rgb[0] * image->hist[i] / scale : 0;
            image->image[3*i + 1] = scale > 0 ? rgb[1] * image->hist[i] / scale : 0;
            image->image[3*i + 2] = scale > 0 ? rgb[2] * image->hist[i] / scale : 0;
#endif
            i++;
        }
    }
}

void capture_massive(dist_t radius, struct massive_particle *P, size_t NP, struct image *image, int clear, int rgb[3])
{
    int32_t r, c;

    if (clear) 
    {
        memset(image->image, 0, 3*image->nc*image->nr*sizeof(*(image->image)));
        memset(image->hist,  0, 1*image->nc*image->nr*sizeof(*(image->hist)));
    }

    size_t i;
    for (i=0; i < NP; i++)
    {
        c = ( P[i].x[0] + radius) / (2.0*radius) * image->nc;
        r = (-P[i].x[1] + radius) / (2.0*radius) * image->nr;
        if (!(0 <= r && r < image->nc)) continue;
        if (!(0 <= c && c < image->nr)) continue;
        image->image[3*(r*image->nc + c) + 0] = rgb[0];
        image->image[3*(r*image->nc + c) + 1] = rgb[1];
        image->image[3*(r*image->nc + c) + 2] = rgb[2];
    }
}


void show_image_formats()
{
    int i;
    for (i=0; supported_image_formats[i].name != NULL; i++)
    {
        eprintf("%s - %s\n", supported_image_formats[i].name, supported_image_formats[i].desc);
    }
}

int check_image_format(char *f)
{
    int i;
    for (i=0; supported_image_formats[i].name != NULL; i++)
    {
        if (!strcmp(f, supported_image_formats[i].name))
            return 1;
    }
    return 0;
}

//==============================================================================
//                               save_image
//==============================================================================
void save_image(int step, struct image *image)
{
    char *fname = make_message("%s%05i.%s", image->name, step, image->format);
    if (!strcmp("png", image->format))
        save_image_png(fname, image->image, image->nr, image->nc);
    else
    {
        // We should already have checked that the user specified format was supported.
        assert(0); 
    }
    free(fname);
}

//==============================================================================
//                               save_image_png
//==============================================================================
int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols)
{
    size_t i;
    FILE *fp = fopen(fname, "wb");

    if (fp == NULL)
    {
        fprintf(stderr, "Can't open %s\n", fname);
        return 1;
    }

    png_structp png_ptr = png_create_write_struct
       (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) return 1;

    png_init_io(png_ptr, fp);

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
       png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
       return 1;
    }

    png_set_IHDR(png_ptr, info_ptr, ncols, nrows,
           8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
           PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);

    uint32_t row_stride = ncols * 3;

    for (i=0; i < nrows; i++)
    {
        png_bytep row_pointer = & (image[i * row_stride]);
        png_write_row(png_ptr, row_pointer);
    }

    png_write_end(png_ptr, info_ptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(fp);

    return 0;
}

void show_snapshot_formats()
{
    int i;
    for (i=0; supported_snapshot_formats[i].name != NULL; i++)
    {
        eprintf("%s - %s\n", supported_snapshot_formats[i].name, 
                             supported_snapshot_formats[i].desc);
    }
}

int check_snapshot_format(char *f)
{
    int i;
    for (i=0; supported_snapshot_formats[i].name != NULL; i++)
    {
        if (!strcmp(f, supported_snapshot_formats[i].name))
            return 1;
    }
    return 0;
}

int save_snapshot_ascii(char *fname, struct binary_header *hdr, struct particle *P, size_t NP, char *col_sep)
{
    size_t i;
    int ret_code = 1;

    FILE *fp;
    
    if (!strcmp(fname, "-"))
        fp = stdout;
    else
        fp = fopen(fname, "wt");

    if (!fp)
    {
        errmsg("Can't open ascii file for writing (%s)", fname);
        ret_code = 0;
        goto cleanup;
    }

    if (fp != stdout)
    {
        time_t *t = (time_t *)&(hdr->creation_time);
        fprintf(fp, 
        "# ASCII snapshot\n"
        "# " GRANDVAL_FULL_PROGRAM_NAME "\n"
        "# Generated on %s"
        "# x y z vx vy vz E id\n",
        asctime(localtime(t)));
    }

    for (i=0; i < NP; i++)
    {
#define F "%24.15e"
        fprintf(fp, F" "F" "F" "F" "F" "F"\n", 
            P[i].x[0], P[i].x[1], P[i].x[2],
            P[i].v[0], P[i].v[1], P[i].v[2], P[i].energy_pp, P[i].id);
    }

    if (fp != stdout)
        fclose(fp);

cleanup:

    return ret_code;
}

int save_snapshot_binary(char *fname, struct binary_header *hdr, struct particle *P, size_t NP, char *col_sep)
{
    size_t i;
    int ret_code = 1;
    FILE *fp;
    
    if (!strcmp(fname, "-"))
        fp = stdout;
    else
        fp = fopen(fname, "wb");

    if (!fp)
    {
        errmsg("Can't open binary file for writing (%s)", fname);
        ret_code = 0;
        goto cleanup;
    }

    fwrite(hdr, sizeof(*hdr), 1, fp);
    fseek(fp, hdr->sizeof_header, SEEK_SET);

    for (i=0; i < NP; i++)
    {
        fwrite(P[i].x, sizeof(*P[i].x), 3, fp);
        fwrite(P[i].v, sizeof(*P[i].v), 3, fp);
        fwrite(&P[i].energy_pp, sizeof(P[i].energy_pp), 1, fp);
        fwrite(&P[i].id, sizeof(P[i].id), 1, fp);
        
    }

    if (fp != stdout)
        fclose(fp);

cleanup:
    return ret_code;
}

int save_snapshot(int step, struct io *io, struct binary_header *hdr, struct particle *P, size_t NP)
{
    int ret_code = 1;
    char *fname;
    
    if (!strcmp(io->name, "-"))
        fname = make_message("-");
    else
    {
        fname = make_message("%s%05i.%s", io->name, step, io->format);

        FILE *fp = fopen(fname, "r+");
        if (fp != NULL && !io->overwrite)
        {
            errmsg("File already exists and overwriting was not allowed (%s).", fname);
            ret_code = 0;
            goto cleanup;
        }
        if (fp != NULL)
            fclose(fp);
        eprintf("Writing to %s\n", fname);
    }

    if (!strcmp("ascii", io->format))
        ret_code = save_snapshot_ascii(fname, hdr, P, NP, " ");
    else if (!strcmp("binary", io->format))
        ret_code = save_snapshot_binary(fname, hdr, P, NP, " ");
    else
        assert(0);

cleanup:

    free(fname);

    return ret_code;
}

void show_binary_snapshot(char *fname)
{
    size_t i;
    FILE *fp;

    if (!strcmp(fname, "-"))
        fp = stdin;
    else
        fp = fopen(fname, "rb");

    if (!fp)
    {
        errmsg("Can't open binary file for writing (%s)", fname);
        return;
    }

    struct binary_header hdr;
    struct particle P;

    if (fread(&hdr, sizeof(char), sizeof(hdr), fp) < sizeof(hdr))
    {
        errmsg("Unexpected end of binary file while reading header.");
        return;
    }

    hdr.major_version = le16toh(hdr.major_version);
    hdr.minor_version = le16toh(hdr.minor_version);
    hdr.sizeof_header = le32toh(hdr.sizeof_header);
    //hdr.sizeof_pos_t  = (uint8_t)sizeof(pos_t);
    //hdr.sizeof_vel_t  = (uint8_t)sizeof(vel_t);
    hdr.creation_time   = le64toh(hdr.creation_time);
    hdr.simulation_step = le64toh(hdr.simulation_step);
    //hdr.simulation_time = sim_time;
    hdr.output_index    = le64toh(hdr.output_index);

    if (fseek(fp, hdr.sizeof_header, SEEK_SET) == -1)
    {
        errmsg("Unable to move to particle data in file. The file is possibly too short.");
        errmsg(strerror(errno));
        return;
    }

    time_t *t = (time_t *)&(hdr.creation_time);

    printf("# Binary file %s\n", fname);
    printf("# grandval v%i.%i\n", hdr.major_version, hdr.minor_version);
    printf("# Generated on %s", asctime(localtime(t)));
    printf("# Step %ld\n", hdr.simulation_step);
    printf("# Time %f\n", hdr.simulation_time);
    printf("# Output %ld\n", hdr.output_index);
    printf("# %ld particles\n", hdr.Nparticles);
    printf("# This file format:\n");
    printf("# x y z vx vy vz E id\n");

#define READN(fp, type, n, lhs) do { \
    size_t _i; \
    type _x[(n)]; \
    fread(_x, sizeof(type), (n), (fp)); \
    for (_i=0; _i < (n); _i++) (lhs)[_i] = _x[_i]; \
} while (0)


    for (i=0; i < hdr.Nparticles; i++)
    {
        if (hdr.sizeof_pos_t == sizeof(float))
            READN(fp, float, 3, P.x);
        else if (hdr.sizeof_pos_t == sizeof(double))
            READN(fp, double, 3, P.x);
        else
            assert(0);

        if (hdr.sizeof_vel_t == sizeof(float))
            READN(fp, float, 3, P.v);
        else if (hdr.sizeof_vel_t == sizeof(double))
            READN(fp, double, 3, P.v);
        else
            assert(0);

        if (hdr.sizeof_energy_t == sizeof(float))
            READN(fp, float, 1, &P.energy_pp);
        else if (hdr.sizeof_energy_t == sizeof(double))
            READN(fp, double, 1, &P.energy_pp);
        else
            assert(0);
            
        READN(fp, int, 1, &P.id);



#define F "%24.15e"
        printf(F" "F" "F" "F" "F" "F" "F" %i\n", 
            P.x[0], P.x[1], P.x[2],
            P.v[0], P.v[1], P.v[2], P.energy_pp, P.id);
    }


    if (fp != stdin)
        fclose(fp);
}

