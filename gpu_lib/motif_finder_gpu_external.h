#ifndef MOTIF_FINDER_MOTIFFINDERGPUEXTERNAL_H
#define MOTIF_FINDER_MOTIFFINDERGPUEXTERNAL_H

#include <vector>
#include <stdint.h>


struct GpuExternalParams {
    GpuExternalParams(const std::vector<uint32_t> &motif_hashes)
        : motif_hashes(motif_hashes)
    {}

    std::vector<std::string> sequences;
    bool search_complementary;
    const std::vector<uint32_t> &motif_hashes;
    uint32_t gpu_count { 0 };
    bool unified_memory {UNIFIED_MEMORY_SUPPORT};
    uint32_t motif_range_size {MOT_PER_CHUNK};
    uint32_t threads_per_block {THREADS_PER_BLOCK};
};

void find_motifs_external_gpu(std::vector<uint16_t> &out_motif_weights, const GpuExternalParams &params);

#endif //MOTIF_FINDER_MOTIFFINDERGPUEXTERNAL_H
