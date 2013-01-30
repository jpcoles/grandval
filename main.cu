#include <cuda.h>
#include <getopt.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "grandval.h"
#include "nbody.h"
#include "ic.h"
#include "io.h"

#define WITH_CUDA 1

size_t cuda_threads_per_block;
size_t cuda_blocks_x, cuda_blocks_y;


#if 0
void step(struct potential *phi, struct particle *p, tyme_t dt)
{
    int i;
    acc_t a[3];

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    phi->accel(phi, p, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}
#endif

void write_positions(struct particle *P, size_t NP, int first_output)
{
    size_t i;
    char *mode;
    char fname[256];

    if (first_output)
        mode = "w";
    else
        mode = "a+";

    FILE *fp = stdout;

    for (i=0; i < NP; i++)
    {
        //sprintf(fname, "gv-%05i.out", i);
        //FILE *fp = fopen(fname, mode);
#if WITH_INTEGERS
        fprintf(fp, "%ld %ld %ld\n", P[i].x[0], P[i].x[1], P[i].x[2]);
#else
        fprintf(fp, "%f %f %f\n", P[i].x[0], P[i].x[1], P[i].x[2]);
#endif
        //fclose(fp);
    }
}

#if 0
void step_all(struct potential *phi, struct particle *P, size_t NP, tyme_t dt)
{
    size_t i;
    #pragma omp parallel for
    for (i=0; i < NP; i++)
        step(phi, P + i, dt);
}
#endif

void show_cuda_devices()
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    int device;
    if (deviceCount == 0)
    {
        printf("No CUDA devices found.\n");
    }

    for (device = 0; device < deviceCount; ++device) 
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        printf("Device %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
        printf("              totalGlobalMem %ld (approx. %ld particles).\n", deviceProp.totalGlobalMem, deviceProp.totalGlobalMem / (sizeof(struct particle)));
        printf("              maxThreadsPerMultiProcessor %d.\n", deviceProp.maxThreadsPerMultiProcessor);
        printf("              maxThreadsPerBlock %d.\n", deviceProp.maxThreadsPerBlock);
        printf("              maxThreadsDim %d,%d,%d.\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
        printf("              maxGridSize %d,%d,%d.\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
        printf("              warpSize %d\n", deviceProp.warpSize);
        printf("              regsPerBlock %d\n", deviceProp.regsPerBlock);

        cuda_threads_per_block = deviceProp.maxThreadsPerBlock;
        cuda_threads_per_block = 512;
        cuda_blocks_x = deviceProp.maxGridSize[0];
        cuda_blocks_y = deviceProp.maxGridSize[1];
    }
}

int select_cuda_device(int device)
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
    {
        printf("No CUDA devices found.\n");
        return 0;
    }

    if (device >= deviceCount)
        return 0;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    cudaError_t err = cudaGetLastError();
    if (err)
    {
        eprintf("Kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    printf("Device %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
    printf("              totalGlobalMem %ld (approx. %ld particles).\n", deviceProp.totalGlobalMem, deviceProp.totalGlobalMem / (sizeof(struct particle)));
    printf("              maxThreadsPerMultiProcessor %d.\n", deviceProp.maxThreadsPerMultiProcessor);
    printf("              maxThreadsPerBlock %d.\n", deviceProp.maxThreadsPerBlock);
    printf("              maxThreadsDim %d,%d,%d.\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
    printf("              maxGridSize %d,%d,%d.\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
    printf("              warpSize %d\n", deviceProp.warpSize);
    printf("              regsPerBlock %d\n", deviceProp.regsPerBlock);

    cuda_threads_per_block = deviceProp.maxThreadsPerBlock;
    cuda_threads_per_block = 512;
    cuda_blocks_x = deviceProp.maxGridSize[0];
    cuda_blocks_y = deviceProp.maxGridSize[1];

    return 1;
}

void show_devices()
{
    show_cuda_devices();
}

void usage()
{
    eprintf("grandval\n");
}

void parse_command_line(int argc, char **argv, struct program_options *opt)
{
    int option_index = 0;
    static struct option long_options[] = 
    {
        {"seed",            1, 0, 0},
        {"verbosity",       1, 0, 0},
        {"dt",              1, 0, 0},
        {"ic",              1, 0, 0},
        {"Nimages",         1, 0, 0},
        {"Rimages",         1, 0, 0},
        {"show-ics",        0, 0, 0},
        {"show-devices",    0, 0, 0},
        {"cuda-device",     0, 0, 0},
        {0, 0, 0, 0}
    };

    while (1)
    {
        int c = getopt_long(argc, argv, "N:vR:", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) 
        {
            case 0:
                     if OPTSTR("dt")              SET_OPTION(opt->dt, atof(optarg));
                else if OPTSTR("Tmax")            SET_OPTION(opt->Tmax, atof(optarg));
                else if OPTSTR("Nimages")         SET_OPTION(opt->Nimages, atoi(optarg));
                else if OPTSTR("Rimages")         SET_OPTION(opt->Rimages, atof(optarg));
                else if OPTSTR("ic")              SET_OPTION(opt->ic_name, optarg);

                else if OPTSTR("show-ics")
                {
                    show_initial_conditions();
                    exit(EXIT_SUCCESS);
                }
                else if OPTSTR("show-devices")
                {
                    show_devices();
                    exit(EXIT_SUCCESS);
                }
                else if OPTSTR("cuda-device")
                {
                    SET_OPTION(opt->cuda_device, atoi(optarg));
                    exit(EXIT_SUCCESS);
                }
                break;

            case 'N': SET_OPTION(opt->Nparticles, atol(optarg)); break;
            case 'R': SET_OPTION(opt->R,          atof(optarg)); break;

            case ':':
            case '?':
                exit(EXIT_FAILURE);
            default:
                usage();
                exit(EXIT_FAILURE);
                break;
        }
    }
}


int main(int argc, char **argv)
{
    size_t NP;
    tyme_t Tmax;
    tyme_t t = 0;
    tyme_t dt;
    dist_t R;

    int Ncaptures;

    struct program_options default_opts;
    memset(&default_opts, 0, sizeof(default_opts));


    SET_OPTION(default_opts.Nparticles, 1000);
    SET_OPTION(default_opts.dt,         0.1);
    SET_OPTION(default_opts.Tmax,       50);
    SET_OPTION(default_opts.R,          0);
    SET_OPTION(default_opts.Nimages,    500);
    SET_OPTION(default_opts.ic_name,    "line");

    struct program_options opts;
    memset(&opts, 0, sizeof(opts));

    int curr_step;

    struct nbody *nbody;
    ic_function ic_f;

    int cuda_device=0;

    int red[3] = {255,0,0};
    int grey[3] = {255,255,255};


    struct image image;
    image.nc = 512;
    image.nr = 512;
    image.image  = (unsigned char *)calloc(3 * image.nr * image.nc, sizeof(*image.image));
    image.hist = (int *)calloc(1 * image.nr * image.nc, sizeof(*image.hist));
    dist_t Rcapture;

    parse_command_line(argc, argv, &opts);

    if (opts.ic_name_set)
    {
        if (!find_ic(opts.ic_name, &ic_f))
        {
            eprintf("Unknown initial condition (%s)\n", optarg);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        assert(find_ic(default_opts.ic_name, &ic_f));
    }

#define ASSIGN_OPTION(v, k, o, d) do { if (o.k ## _set) { v=o.k; } else {v = d.k;} } while (0)
    ASSIGN_OPTION(NP, Nparticles, opts, default_opts);
    ASSIGN_OPTION(dt, dt, opts, default_opts);
    ASSIGN_OPTION(cuda_device, cuda_device, opts, default_opts);
    ASSIGN_OPTION(Tmax, Tmax, opts, default_opts);
    ASSIGN_OPTION(Rcapture, Rimages, opts, default_opts);
    ASSIGN_OPTION(Ncaptures, Nimages, opts, default_opts);
    ASSIGN_OPTION(R, R, opts, default_opts);

    eprintf("NP %ld\n", NP);
    eprintf("dt %f\n", dt);
    eprintf("cuda device %i\n", cuda_device);
    eprintf("Tmax %f\n", Tmax);
    eprintf("Rcapture %f\n", Rcapture);
    eprintf("Ncaptures %i\n", Ncaptures);
    eprintf("R %f\n", R);

    if (!select_cuda_device(cuda_device))
        exit(EXIT_FAILURE);

    struct particle *P = (struct particle *)malloc(NP * sizeof(*P));
    assert(P != NULL);

    ic_f(P, NP, R);

    struct potential phi;

    void *phi_data = nbody_init(&phi);
    nbody = (struct nbody *)phi_data;
    nbody_create_potential(nbody, 2);

    phi.set_particles(phi_data, P, NP);

    double t_next_capture = Tmax / Ncaptures;
    int curr_capture = 0;

    phi.step_particles(phi_data, dt);
    phi.advance(phi_data, 0);

    P = 0;
    NP = 0;
    phi.get_particles(phi_data, &P, &NP);

    if (Ncaptures)
    {
        capture(Rcapture, P, NP, &image, 1, grey);
        capture_massive(Rcapture, nbody->phi.P, nbody->phi.N, &image, 0, red);
        save_snapshot(curr_capture, &image);
        curr_capture++;
    }

    for (curr_step = 0, t = dt;
         t < Tmax+dt;
         t += dt, curr_step++)
    {
        if (t > Tmax) t = Tmax;

        phi.step_particles(phi_data, dt);
        phi.advance(phi_data, t);
        phi.get_particles(phi_data, &P, &NP);

        fprintf(stderr, "Rcapture %f\n", Rcapture);

        if (Ncaptures && t > t_next_capture)
        {
            //write_positions(P, NP, curr_step == 0);
            capture(Rcapture, P, NP, &image, 1, grey);
            capture_massive(Rcapture, nbody->phi.P, nbody->phi.N, &image, 0, red);
            save_snapshot(curr_capture, &image);
            curr_capture++;
            t_next_capture += Tmax / Ncaptures;
        }
    }


    phi.get_particles(phi_data, &P, &NP);
    //nbody_get_particles(nbody, &P, &NP);
    if (Ncaptures)
    {
        capture(Rcapture, P, NP, &image, 1, grey);
        capture_massive(Rcapture, nbody->phi.P, nbody->phi.N, &image, 0, red);
        save_snapshot(curr_capture, &image);
    }
    //write_positions(P, NP, curr_step == 0);

    //nbody_free(nbody);
    phi.free(phi_data);

    return 0;
}

