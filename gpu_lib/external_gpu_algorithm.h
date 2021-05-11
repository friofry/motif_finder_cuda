#endif //MOTIF_FINDER_EXTERNAL_CPU_ALGORITHM_H

#ifndef MOTIF_FINDER_MOTIF_FINDER_GPU_EXTERNAL_H
#define MOTIF_FINDER_MOTIF_FINDER_GPU_EXTERNAL_H

#include <vector>

#include "gpu_cuda_params.h"
#include <config.h>
#include <sequence_hashes.h>

// Get motif occurrences and store to out_motif_weights
// @param motif_hashes - motif hashes
// @param sequence_hashes - sequences hashes
// @param out_motif_weights - [output] array for motif occurrences
// @param params - CUDA gpu parameters
void external_gpu_algorithm(const std::vector<uint32_t> &motif_hashes,
                            const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights,
                            const GpuCudaParams &params);

#endif // MOTIF_FINDER_MOTIF_FINDER_GPU_EXTERNAL_H
