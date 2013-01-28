#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include "grandval.h"

#include "io.h"
#include "color_ramp.h"


//==============================================================================
//                                  capture
//==============================================================================
void capture(double radius, struct particle *P, size_t NP, struct image *image, int clear, int rgb[3])
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

void capture_massive(double radius, struct massive_particle *P, size_t NP, struct image *image, int clear, int rgb[3])
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

static char *make_message(const char *fmt, ...)
{
    /* Guess we need no more than 100 bytes. */
    int n, size = 100;
    char *p, *np;
    va_list ap;

    if ((p = (char *)malloc(size)) == NULL)
        return NULL;

    while (1) 
    {
        /* Try to print in the allocated space. */
        va_start(ap, fmt);
        n = vsnprintf(p, size, fmt, ap);
        va_end(ap);
        /* If that worked, return the string. */
        if (n > -1 && n < size)
            return p;
        /* Else try again with more space. */
        if (n > -1)    /* glibc 2.1 */
            size = n+1; /* precisely what is needed */
        else           /* glibc 2.0 */
            size *= 2;  /* twice the old size */
        if ((np = (char *)realloc (p, size)) == NULL) 
        {
            free(p);
            return NULL;
        } else {
            p = np;
        }
    }
}

//==============================================================================
//                               save_snapshot
//==============================================================================
void save_snapshot(int step, struct image *image)
{
    char *fname = make_message("gv-%05i.png", step);
    save_image_png(fname, image->image, image->nr, image->nc);
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
