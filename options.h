#ifndef OPTIONS_H
#define OPTIONS_H

#include "grandval.h"

void usage();
void parse_command_line(int argc, char **argv, struct program_options *opt);
void show_options(struct program_options *opts, struct program_options *default_opts);
char *make_binary_options_text(struct program_options *opts, struct program_options *default_opts);

#endif
