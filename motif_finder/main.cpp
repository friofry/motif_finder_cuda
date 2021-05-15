#include <iostream>
#include <cstdio>

#include <argo_cuda_params_reader.h>
#include <find_motifs_iterative.h>
#include <fst_reader.h>

#include "create_stat_model.h"

using namespace std;

const char *_default_ini_filename = "init.ini";


int main(int argc, char** argv)
{
    std::ignore = argc;
    std::ignore = argv;
    ArgoCudaParams params = read_ini_file(_default_ini_filename);
    print_argo_cuda_params(params);

    auto sequences = read_fasta(params.positive_sequences.c_str());

    auto stat_model = create_stat_model(params.use_real_nucl_frequences, params.markov_level, sequences, params.complementary);


    find_motifs_iterative();

    return 0;
}
