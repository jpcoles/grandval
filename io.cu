#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "grandval.h"

#include "io.h"
#include "util.h"
#include "color_ramp.h"

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
    {NULL, NULL} 
};

//==============================================================================
//                                  capture
//==============================================================================
void capture(dist_t radius, struct particle *P, size_t NP, struct image *image, int clear, int rgb[3])
{
    int32_t r, c;
    double scale = 0;

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
        if (image->hist[r*image->nc + c] > scale)
            scale = image->hist[r*image->nc + c];
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

int save_snapshot_ascii(char *fname, struct particle *P, size_t NP, char *col_sep)
{
    size_t i;
    int ret_code = 1;

    FILE *fp = fopen(fname, "wt");
    if (!fp)
    {
        errmsg("Can't open file for writing (%s)", fname);
        ret_code = 0;
        goto cleanup;
    }

    for (i=0; i < NP; i++)
    {
        fprintf(fp, "x\n");
    }

    fclose(fp);

cleanup:

    return ret_code;
}

int save_snapshot(int step, struct io *io, struct particle *P, size_t NP)
{
    int ret_code = 1;
    char *fname = make_message("%s%05i.%s", io->name, step, io->format);

    FILE *fp = fopen(fname, "r+");
    if (fp != NULL && io->overwrite)
    {
        errmsg("File already exists and overwriting is was not allowed (%s).", fname);
        ret_code = 0;
        goto cleanup;
    }
    fclose(fp);

    eprintf("Writing to %s\n", fname);

    if (!strcmp("ascii", io->format))
        return save_snapshot_ascii(fname, P, NP, " ");
    if (!strcmp("csv", io->format))
        return save_snapshot_ascii(fname, P, NP, ",");
    else
        assert(0);

cleanup:

    free(fname);

    return ret_code;
}
