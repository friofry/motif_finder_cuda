#include <cstdio>
#include <iostream>

// common_lib
#include <fst_reader.h>

// cpu_lib
#include <external_gpu_algorithm.h>

// iterative_finder_lib
#include <argo_cuda_params_reader.h>
#include <iterative_finder.h>
#include <external_gpu_algorithm.h>

using namespace std;

const char *_default_ini_filename = "init.ini";

int main(int, char **)
{
    // 1. Read parameters
    ArgoCudaParams params = read_ini_file(_default_ini_filename);

    // 2. Prepare algorithm
    GpuCudaParams gpu_params;
    auto cpu_external_algorithm = [&gpu_params](const std::vector<uint32_t> &motif_hashes,
                                     const SequenceHashes &sequence_hashes,
                                     std::vector<uint16_t> &out_motif_weights) {
        external_gpu_algorithm(motif_hashes, sequence_hashes, out_motif_weights, gpu_params);
    };

    // 3. Run iterative search
    find_important_motifs(params, cpu_external_algorithm);

    return 0;
}
