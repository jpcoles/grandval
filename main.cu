#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "grandval.h"
#include "ic.h"
#include "io.h"
#include "options.h"
#include "devices.h"
#include "util.h"
#include "rand.h"

#include "nbody.h"
#include "bar.h"
#include "plummer.h"
#include "hernquist.h"
#include "jaffe.h"
#include "plummer_lumps.h"

size_t phase_space_number_density(struct particle *P, size_t N, dist_t R)
{
    pos_t x[6];
    x[0] = P[0].x[0];
    x[1] = P[0].x[1];
    x[2] = P[0].x[2];
    x[3] = P[0].v[0];
    x[4] = P[0].v[1];
    x[5] = P[0].v[2];
    const dist_t R2 = R*R;

    size_t i,j;
    size_t n=0;
    for (i=0; i < N; i++)
    {
        dist_t d[6];
        for (j=0; j < 3; j++)
        {
            d[j+0] = P[i].x[j] - x[j+0];
            d[j+3] = P[i].v[j] - x[j+3];
        }
        dist_t r2 = 0;
        for (j=0; j < 6; j++)
            r2 += d[j] * d[j];

        n += r2 <= R2;
        //if (r2 > R2) eprintf("%24.15f %24.15f\n", r2, R2);
    }

    return n;
}

int main(int argc, char **argv)
{
    size_t NP;
    tyme_t Tmax;
    tyme_t t = 0;
    tyme_t dt;
    dist_t R;

    dist_t psR = 1e-6;

    int ret_code = EXIT_SUCCESS;

    FILE *energy_fp = NULL;
    struct particle *P = NULL;
    void *phi_data = NULL;

    struct potential phi;
    struct program_options default_opts;
    struct program_options opts;

    struct binary_header *hdr = NULL;

    int curr_step;

    struct nbody *nbody;
    ic_function ic_f;

    int red[3] = {255,0,0};
    int grey[3] = {255,255,255};


    struct image image;
    struct io io;

    double t_next_capture=0;
    int curr_capture=0;

    double t_next_save=0;
    int curr_save=0;

    memset(&phi, 0, sizeof(phi));
    memset(&default_opts, 0, sizeof(default_opts));
    memset(&opts, 0, sizeof(opts));

    atexit(print_errmsg);

    SET_OPTION(default_opts.Nparticles,             1000);
    SET_OPTION(default_opts.dt,                     0.1);
    SET_OPTION(default_opts.Tmax,                   50);
    SET_OPTION(default_opts.R,                      10);
    SET_OPTION(default_opts.Nimages,                500);
    SET_OPTION(default_opts.Rimages,                10);
    SET_OPTION(default_opts.ic_name,                "line");
    SET_OPTION(default_opts.potential_name,         "nbody");
    SET_OPTION(default_opts.Nsnapshots,             0);
    SET_OPTION(default_opts.snapshot_name,          "gv-");
    SET_OPTION(default_opts.snapshot_format,        "ascii");
    SET_OPTION(default_opts.image_name,             "gv-");
    SET_OPTION(default_opts.image_format,           "png");
    SET_OPTION(default_opts.overwrite,              0);
    SET_OPTION(default_opts.random_seed,            time(NULL));

    nbody_init(&phi);       add_potential(&phi);
    bar_init(&phi);         add_potential(&phi);
    plummer_init(&phi);     add_potential(&phi);
    hernquist_init(&phi);   add_potential(&phi);
    jaffe_init(&phi);       add_potential(&phi);
    plummer_lumps_init(&phi); add_potential(&phi);

    parse_command_line(argc, argv, &opts);

#define ASSIGN_OPTION(v, k, o, d) do { if (o.k ## _set) { v=o.k; } else {v = d.k;} } while (0)
    ASSIGN_OPTION(NP, Nparticles, opts, default_opts);
    ASSIGN_OPTION(dt, dt, opts, default_opts);
    ASSIGN_OPTION(Tmax, Tmax, opts, default_opts);
    ASSIGN_OPTION(R, R, opts, default_opts);

    ASSIGN_OPTION(opts.random_seed, random_seed, opts, default_opts);

    rand_seed(opts.random_seed);

    if (opts.Nimages_set)
    {
        if (!opts.Rimages_set)
        {
            errmsg("Option --Rimages must be set to save images.");
            goto fail;
        }

        assert(opts.Rimages_set);
        ASSIGN_OPTION(image.name,   image_name,   opts, default_opts);
        ASSIGN_OPTION(image.format, image_format, opts, default_opts);
        if (!check_image_format(image.format))
        {
            errmsg("Image format not supported (%s).", image.format);
            goto fail;
        }

        image.nc = 512;
        image.nr = 512;
        image.image  = (unsigned char *)calloc(3 * image.nr * image.nc, sizeof(*image.image));
        image.hist = (int *)calloc(1 * image.nr * image.nc, sizeof(*image.hist));

        t_next_capture = Tmax / opts.Nimages;
        curr_capture = 0;
    }

    if (opts.Nsnapshots_set)
    {
        if (opts.Nsnapshots < 0)
        {
            errmsg("Number of snapshots requested must not be negative (%ld).", opts.Nsnapshots);
            goto fail;
        }

        ASSIGN_OPTION(io.name,   snapshot_name,   opts, default_opts);
        ASSIGN_OPTION(io.format, snapshot_format, opts, default_opts);
        ASSIGN_OPTION(io.overwrite, overwrite, opts, default_opts);
        if (!check_snapshot_format(io.format))
        {
            errmsg("Snapshot format not supported (%s).", io.format);
            goto fail;
        }

        t_next_save = Tmax / opts.Nsnapshots;
        curr_save = 0;
    }

    if (opts.ic_name_set)
    {
        if (!find_ic(opts.ic_name, &ic_f))
        {
            errmsg("Unknown initial condition (%s).", opts.ic_name);
            goto fail;
        }
    }
    else
    {
        assert(find_ic(default_opts.ic_name, &ic_f));
    }

    if (opts.potential_name_set)
    {
        if (!find_potential(opts.potential_name, &phi))
        {
            errmsg("Unknown potential (%s).", opts.potential_name);
            goto fail;
        }
    }
    else
    {
        assert(find_potential(default_opts.potential_name, &phi));
    }

    if (opts.energy_fname_set)
    {
        if (!opts.overwrite)
        {
            energy_fp = fopen(opts.energy_fname, "r");
            if (energy_fp)
            {
                fclose(energy_fp);
                errmsg("File already exists and overwriting was not allowed (%s).", opts.energy_fname);
                goto fail;
            }
        }

        energy_fp = fopen(opts.energy_fname, "wt");
        if (!energy_fp)
        {
            errmsg("Unable to open energy file (%s).", opts.energy_fname);
            goto fail;
        }
    }

    show_options(&opts, &default_opts);

    if (!select_cuda_device(opts.cuda_device))
        goto fail;

    P = (struct particle *)malloc(NP * sizeof(*P));
    assert(P != NULL);

    ic_f(P, NP, R);

    //printf("PSD %ld\n", phase_space_number_density(P, NP, psR));

    phi_data = phi.create(2);
    nbody = (struct nbody *)phi_data;

    phi.set_particles(phi_data, P, NP);

    if (energy_fp && phi.energy != NULL)
        fprintf(energy_fp, "%24.15e\n", phi.energy(phi_data, P, NP));

    P = 0;
    NP = 0;
    phi.get_particles(phi_data, &P, &NP);

    if (opts.Nimages_set && opts.Nimages)
    {
        capture(opts.Rimages, P, NP, &image, 1, grey);
        //capture_massive(opts.Rimages, nbody->phi.P, nbody->phi.N, &image, 0, red);
        save_image(curr_capture, &image);
        curr_capture++;
    }

    if (opts.Nsnapshots_set && opts.Nsnapshots > 0)
    {
        create_io_header(&hdr, t, curr_save, curr_step, &opts, &default_opts, NP);
        if (!save_snapshot(curr_save, &io, hdr, P, NP))
            goto fail;
        curr_save++;
    }

    curr_step = 0;

    while (1)
    {
        curr_step++;
        t = dt*curr_step;
        if (t > Tmax+dt) break;

        if (t > Tmax) t = Tmax;

        eprintf("Simulation time %f\n", t);

        if (!phi.step_particles(phi_data, dt)) goto fail;
        if (!phi.advance(phi_data, t)) goto fail;
        if (!phi.get_particles(phi_data, &P, &NP)) goto fail;

        if (opts.Nimages_set && opts.Nimages > 0 && t > t_next_capture)
        {
            //write_positions(P, NP, curr_step == 0);
            capture(opts.Rimages, P, NP, &image, 1, grey);
            //capture_massive(opts.Rimages, nbody->phi.P, nbody->phi.N, &image, 0, red);
            save_image(curr_capture, &image);
            curr_capture++;
            t_next_capture += Tmax / opts.Nimages;
        }

        if (opts.Nsnapshots_set && opts.Nsnapshots > 0 && t > t_next_save)
        {
            create_io_header(&hdr, t, curr_save, curr_step, &opts, &default_opts, NP);
            if (!save_snapshot(curr_save, &io, hdr, P, NP))
                goto fail;

            curr_save++;
            t_next_save += Tmax / opts.Nsnapshots;
        }

        //printf("PSD %ld\n", phase_space_number_density(P, NP, psR));

        if (energy_fp && phi.energy != NULL)
            fprintf(energy_fp, "%24.15e\n", phi.energy(phi_data, P, NP));
    }


    phi.get_particles(phi_data, &P, &NP);
    if (opts.Nimages_set && opts.Nimages > 0 && t > t_next_capture)
    {
        capture(opts.Rimages, P, NP, &image, 1, grey);
        //capture_massive(opts.Rimages, nbody->phi.P, nbody->phi.N, &image, 0, red);
        save_image(curr_capture, &image);
    }

    if (opts.Nsnapshots_set && opts.Nsnapshots > 0 && t > t_next_save)
    {
        create_io_header(&hdr, t, curr_save, curr_step, &opts, &default_opts, NP);
        if (!save_snapshot(curr_save, &io, hdr, P, NP))
            goto fail;
    }

    goto cleanup;

fail:
    //print_errmsg();
    ret_code = EXIT_FAILURE;

cleanup:

    if (phi.free) phi.free(phi_data);

    return ret_code;
}

