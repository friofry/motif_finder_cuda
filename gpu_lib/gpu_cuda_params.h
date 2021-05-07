#ifndef MOTIF_FINDER_GPU_CUDA_PARAMS_H
#define MOTIF_FINDER_GPU_CUDA_PARAMS_H
#include <config.h>

struct GpuCudaParams {
    bool unified_memory {UNIFIED_MEMORY_SUPPORT};
    uint32_t motif_range_size {MOT_PER_CHUNK};
    uint32_t threads_per_block {THREADS_PER_BLOCK};
};

#endif //MOTIF_FINDER_GPU_CUDA_PARAMS_H
