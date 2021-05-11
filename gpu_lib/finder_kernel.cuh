#ifndef MOTIF_FINDER_FINDER_KERNEL_CUH
#define MOTIF_FINDER_FINDER_KERNEL_CUH
#include <sequence_hashes.h>

#include "gpu_memory.h"
#include "gpu_cuda_params.h"


// these are wrapper functions for cuda kernels

// EXTERNAL
class GpuExternalMemory {
public:
    GpuExternalMemory(const GpuCudaParams &params, const SequenceHashes &sequence_hashes);
    ~GpuExternalMemory();

    GpuMemoryPtr alloc;
    uint32_t shared_memory_size {0};
    uint32_t sequences_count {0};
    uint32_t weights_count {0};

    uint16_t *weights_out {nullptr};
    uint32_t *seq_hashes {nullptr};
    uint32_t *seq_lengths {nullptr};
    uint32_t *seq_begins {nullptr};
    uint32_t *motif_hashes {nullptr};
};

void motif_finder_gpu_external(
    const std::vector<uint32_t> &_motif_hashes,
    const GpuExternalMemory &gpu_mem,
    const GpuCudaParams &params,
    std::vector<uint16_t> &out_motif_weights,
    uint32_t motif_idx_offset,
    uint32_t motifs_count,
    int device_id);

// INTERNAL
class GpuInternalMemory {
public:
    GpuInternalMemory(const GpuCudaParams &params, const SequenceHashes &sequence_hashes);
    ~GpuInternalMemory();

    GpuMemoryPtr alloc;
    uint32_t shared_memory_size {0};
    uint32_t sequences_count {0};
    uint32_t weights_count {0};
    uint16_t *weights_out {nullptr};
    uint32_t *seq_hashes {nullptr};
    uint32_t *seq_lengths {nullptr};
    uint32_t *seq_begins {nullptr};
};

void motif_finder_gpu_internal(
    const GpuInternalMemory &gpu_mem,
    const GpuCudaParams &params,
    std::vector<uint16_t> &out_motif_weights,
    uint32_t motif_idx_offset,
    uint32_t motifs_count,
    int device_id);

#endif //MOTIF_FINDER_FINDER_KERNEL_CUH
