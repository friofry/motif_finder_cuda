#include "external_gpu_algorithm.h"

#include <thread>

#include <config.h>
#include <hash_conversions.h>
#include <run_parallel.h>
#include <safe_counter.h>

#include <gpu_info.h>

#include "finder_kernel.cuh"

void external_gpu_algorithm(const std::vector<uint32_t> &motif_hashes,
                            const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights,
                            const GpuCudaParams &params)
{
    out_motif_weights.resize(motif_hashes.size(), 0);
    uint32_t threads = (params.gpu_count > 0) ? params.gpu_count : 1;

    SafeCounter motifs_counter(motif_hashes.size());

    run_parallel(threads, [&](uint32_t thread_id) {
        GpuExternalMemory gpu_memory(params, sequence_hashes);
        while (true) {
            const auto range = motifs_counter.get_and_increment_range_info(params.motif_range_size);
            if (range.count() == 0) {
                break;
            }
            motif_finder_gpu_external(
                motif_hashes, gpu_memory, params, out_motif_weights, range.start, range.count(), thread_id);
        }
    });
}