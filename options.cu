#include <stdio.h>
#include <getopt.h>
#include "grandval.h"
#include "options.h"
#include "potential.h"
#include "ic.h"
#include "devices.h"
#include "io.h"
#include "util.h"

void usage()
{
    eprintf(
    "Usage: grandval OPTIONS\n"
    "\n"
    "where OPTIONS can be any of\n"
    "   --Tmax <float>                      Set total simulation time.\n"
    "   --dt <float>                        Set the time step of the test particles.\n"
    "   --ic <string>                       Choose the initial particle conditions.\n"
    "   --show-ics                          Show the available initial conditions.\n"
    "   --potential <string>                Choose a potential to evolve the particles in.\n"
    "   --show-potentials                   Show the available potentials.\n"
    "   --Nsnapshots <float>                Set the number of snapshots.\n"
    "   --snapshot-name <string>            Set the filename prefix of the snapshot outputs.\n"
    "   --snapshot-format <string>          Set the snapshot format.\n"
    "   --show-snapshot-formats             Show the available snapshot formats.\n"
    "   --overwrite                         Overwrite any existing output files.\n"
    "   --Nimages <int>                     Set the number of image outputs.\n"
    "   --Rimages <float>                   Set the physical radius of the image outputs.\n"
    "   --image-name <string>               Set the filename prefix of the image outputs.\n"
    "   --image-format <string>             Set the image format.\n"
    "   --show-image-formats                Show the available image formats\n"
    "   --cuda-device <int>                 Set the CUDA compatible device.\n"              
    "   --show-devices                      Show the available devices.\n"
    "   --seed <int>                        Set the random number generator seed.\n"
    "   --save-energy <string>              Save the total energy at each time step to the specified file.\n"
    "   --show-binary <string>              Display a binary snapshot in ASCII text and quit.\n"
    "\n"
    "Send bug reports to Jonathan Coles <jonathan@physik.uzh.ch>\n"
    );
}

void parse_command_line(int argc, char **argv, struct program_options *opt)
{
    int option_index = 0;
    static struct option long_options[] = 
    {
        {"help",            0, 0, 0},
        {"seed",            1, 0, 0},
        {"verbosity",       1, 0, 0},
        {"Tmax",            1, 0, 0},
        {"dt",              1, 0, 0},
        {"ic",              1, 0, 0},
        {"show-ics",        0, 0, 0},
        {"potential",       1, 0, 0},
        {"show-potentials",        0, 0, 0},
        {"Nsnapshots",         1, 0, 0},
        {"snapshot-name",         1, 0, 0},
        {"snapshot-format",         1, 0, 0},
        {"show-snapshot-formats",  0, 0, 0},
        {"Nimages",         1, 0, 0},
        {"Rimages",         1, 0, 0},
        {"image-name",     1, 0, 0},
        {"image-format",     1, 0, 0},
        {"show-image-formats",     0, 0, 0},
        {"show-devices",    0, 0, 0},
        {"cuda-device",     1, 0, 0},
        {"overwrite",     0, 0, 0},
        {"save-energy",     1, 0, 0},
        {"show-binary",     1, 0, 0},
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
                     if OPTSTR("help")            { usage(); exit(EXIT_SUCCESS); }
                else if OPTSTR("seed")            SET_OPTION(opt->random_seed,          atol(optarg));
                else if OPTSTR("dt")              SET_OPTION(opt->dt,                   atof(optarg));
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
                else if OPTSTR("overwrite")       SET_OPTION(opt->overwrite,                      1);
                else if OPTSTR("cuda-device")     SET_OPTION(opt->cuda_device,          atoi(optarg));
                else if OPTSTR("save-energy")     SET_OPTION(opt->energy_fname,              optarg);

                else if OPTSTR("show-potentials") { show_potentials();         exit(EXIT_SUCCESS); }
                else if OPTSTR("show-ics")        { show_initial_conditions(); exit(EXIT_SUCCESS); }
                else if OPTSTR("show-devices")    { show_devices();            exit(EXIT_SUCCESS); }
                else if OPTSTR("show-image-formats")    { show_image_formats();            exit(EXIT_SUCCESS); }
                else if OPTSTR("show-snapshot-formats") { show_snapshot_formats();         exit(EXIT_SUCCESS); }
                else if OPTSTR("show-binary")           { show_binary_snapshot(optarg);    exit(EXIT_SUCCESS); }
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
    PRINT_OPT("Nparticles                %ld", Nparticles);
    PRINT_OPT("dt                        %f", dt);
    PRINT_OPT("Tmax                      %f", Tmax);
    PRINT_OPT("Nimages                   %i", Nimages);
    PRINT_OPT("Rimages                   %f", Rimages);
    PRINT_OPT("Image name                %s", image_name);
    PRINT_OPT("Image format              %s", image_format);
    PRINT_OPT("R                         %f", R);
    PRINT_OPT("IC                        %s", ic_name);
    PRINT_OPT("Potential                 %s", potential_name);
    PRINT_OPT("Cuda Device               %i", cuda_device);
    PRINT_OPT("Nsnapshots                %ld", Nsnapshots);
    PRINT_OPT("Snapshot name             %s", snapshot_name);
    PRINT_OPT("Snapshot format           %s", snapshot_format);
    PRINT_OPT("Overwrite                 %i", overwrite);
    PRINT_OPT("Random seed               %ld", random_seed);
    PRINT_OPT("Energy filename           %s", energy_fname);
}

char *make_binary_options_text(struct program_options *opts, struct program_options *default_opts)
{
    char *msg = (char *)malloc(1024 * sizeof(*msg));

#define BIN_OPT(str, o) do { char *m; if (opts -> o ## _set) m=make_message("S " #o "\n", opts -> o); else m=make_message("D " #o "\n", default_opts -> o); strcat(msg, m); free(m); } while (0)
    BIN_OPT("%ld",    Nparticles);
    BIN_OPT("%f",     dt);
    BIN_OPT("%f",     Tmax);
    BIN_OPT("%i",     Nimages);
    BIN_OPT("%f",     Rimages);
    BIN_OPT("%s",     image_name);
    BIN_OPT("%s",     image_format);
    BIN_OPT("%f",     R);
    BIN_OPT("%s",     ic_name);
    BIN_OPT("%s",     potential_name);
    BIN_OPT("%i",     cuda_device);
    BIN_OPT("%ld",    Nsnapshots);
    BIN_OPT("%s",     snapshot_name);
    BIN_OPT("%s",     snapshot_format);
    BIN_OPT("%i",     overwrite);
    BIN_OPT("%ld",    random_seed);
    BIN_OPT("%s",     energy_fname);

    return msg;
}

