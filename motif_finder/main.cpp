#include <iostream>
#include <cstdio>

// common_lib
#include <fst_reader.h>

// cpu_lib
#include <external_cpu_algorithm.h>

// iterative_finder_lib
#include <argo_cuda_params_reader.h>
#include <iterative_finder.h>

#include "create_stat_model.h"

using namespace std;

const char *_default_ini_filename = "init.ini";

int main(int argc, char** argv)
{
    std::ignore = argc;
    std::ignore = argv;
    ArgoCudaParams params = read_ini_file(_default_ini_filename);
    print_argo_cuda_params(params);

    auto cpu_external_algorithm = [](const std::vector<uint32_t> &motif_hashes,
                                 const SequenceHashes &sequence_hashes,
                                 std::vector<uint16_t> &out_motif_weights) {
        external_cpu_algorithm(motif_hashes, sequence_hashes, out_motif_weights);
    };

    ImportantMotifFinder finder(params, cpu_external_algorithm);
    finder.find();

    return 0;
}
