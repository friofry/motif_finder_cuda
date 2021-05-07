#include <sstream>
#include <vector>
#include <cstdint>

#include <config.h>
#include <motif_finder_gpu_internal.h>
#include <motif_finder_gpu_external.h>
#include <fst_reader.h>

using namespace std;

void run()
{
    GpuInternalParams params;
    params.sequences = read_fasta("test.fst");
    params.search_complementary = false;
    params.gpu_count = 1;
    params.unified_memory = false;
    params.motif_range_size = MOT_PER_CHUNK;
    params.threads_per_block = THREADS_PER_BLOCK;

    std::vector<uint16_t> res;
    find_motifs_internal_gpu(res, params);


    std::vector<uint32_t> motif_hashes;
    GpuExternalParams external_params(motif_hashes);
    find_motifs_external_gpu(res, params);
}

int main(int argc, char **argv)
{
    argc = 0;
    argv = 0;
    run();
    return 0;
}
