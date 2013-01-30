#include <stdio.h>
#include <getopt.h>
#include "grandval.h"
#include "options.h"
#include "potential.h"
#include "ic.h"
#include "devices.h"
#include "io.h"

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
        {"potential",       1, 0, 0},
        {"Nimages",         1, 0, 0},
        {"Rimages",         1, 0, 0},
        {"show-ics",        0, 0, 0},
        {"show-devices",    0, 0, 0},
        {"cuda-device",     1, 0, 0},
        {"image-name",     1, 0, 0},
        {"image-format",     1, 0, 0},
        {"show-image-formats",     0, 0, 0},
        {"Nsnapshots",         1, 0, 0},
        {"snapshot-name",         1, 0, 0},
        {"snapshot-format",         1, 0, 0},
        {"show-snapshot-formats",  0, 0, 0},
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
                     if OPTSTR("dt")              SET_OPTION(opt->dt,                   atof(optarg));
                else if OPTSTR("Tmax")            SET_OPTION(opt->Tmax,                 atof(optarg));
                else if OPTSTR("Nimages")         SET_OPTION(opt->Nimages,              atoi(optarg));
                else if OPTSTR("Rimages")         SET_OPTION(opt->Rimages,              atof(optarg));
                else if OPTSTR("image-name")      SET_OPTION(opt->image_name,                optarg);
                else if OPTSTR("image-format")    SET_OPTION(opt->image_format,              optarg);
                else if OPTSTR("ic")              SET_OPTION(opt->ic_name,                   optarg);
                else if OPTSTR("potential")       SET_OPTION(opt->potential_name,            optarg);
                else if OPTSTR("Nsnapshots")      SET_OPTION(opt->Nsnapshots,           atoi(optarg));
                else if OPTSTR("snapshot-name")   SET_OPTION(opt->snapshot_name,             optarg);
                else if OPTSTR("snapshot-format") SET_OPTION(opt->snapshot_format,           optarg);
                else if OPTSTR("cuda-device")     SET_OPTION(opt->cuda_device,          atoi(optarg));

                else if OPTSTR("show-potentials") { show_potentials();         exit(EXIT_SUCCESS); }
                else if OPTSTR("show-ics")        { show_initial_conditions(); exit(EXIT_SUCCESS); }
                else if OPTSTR("show-devices")    { show_devices();            exit(EXIT_SUCCESS); }
                else if OPTSTR("show-image-formats")    { show_image_formats();            exit(EXIT_SUCCESS); }
                else if OPTSTR("show-snapshot-formats") { show_snapshot_formats();            exit(EXIT_SUCCESS); }
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

void show_options(struct program_options *opts, struct program_options *default_opts)
{
#define PRINT_OPT(str, o) do { if (opts -> o ## _set) eprintf(str "\n", opts -> o); else eprintf(str " (default)\n", default_opts -> o); } while (0)
    PRINT_OPT("Nparticles               %ld", Nparticles);
    PRINT_OPT("dt                        %f", dt);
    PRINT_OPT("Tmax                      %f", Tmax);
    PRINT_OPT("Nimages                   %i", Nimages);
    PRINT_OPT("Rimages                   %f", Rimages);
    PRINT_OPT("Image name                %f", image_name);
    PRINT_OPT("Image format              %f", image_format);
    PRINT_OPT("R                         %f", R);
    PRINT_OPT("IC                        %s", ic_name);
    PRINT_OPT("Potential                 %s", potential_name);
    PRINT_OPT("Cuda Device               %i", cuda_device);
    PRINT_OPT("Nsnapshots                %i", Nsnapshots);
    PRINT_OPT("Snapshot name             %s", snapshot_name);
    PRINT_OPT("Snapshot format           %s", snapshot_format);
}

