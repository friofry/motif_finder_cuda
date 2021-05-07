#ifndef MOTIF_FINDER_MOTIF_FINDER_GPU_INTERNAL_H
#define MOTIF_FINDER_MOTIF_FINDER_GPU_INTERNAL_H

#include <string>
#include <vector>
#include <config.h>

struct GpuInternalParams {
    std::vector<std::string> sequences;
    bool search_complementary { false };
    uint32_t gpu_count { 0 };
    bool unified_memory {UNIFIED_MEMORY_SUPPORT};
    uint32_t motif_range_size {MOT_PER_CHUNK};
    uint32_t threads_per_block {THREADS_PER_BLOCK};
};

void find_motifs_internal_gpu(std::vector<uint16_t> &out_motif_weights, const GpuInternalParams &params);

#endif // MOTIF_FINDER_MOTIF_FINDER_GPU_INTERNAL_H
