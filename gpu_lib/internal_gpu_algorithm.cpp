#include "internal_gpu_algorithm.h"

#include <config.h>
#include <hash_conversions.h>
#include <run_parallel.h>
#include <safe_counter.h>

#include <gpu_info.h>

#include "finder_kernel.cuh"

void internal_gpu_algorithm(const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights,
                            const GpuCudaParams &params)
{
    out_motif_weights.resize(TOTAL_MOT, 0);
    threads = (threads > 0) ? threads : std::thread::hardware_concurrency();

    SafeCounter motifs_counter(TOTAL_MOT);

    run_parallel(threads, [&](uint32_t thread_id) {
        GpuInternalMemory gpu_memory(params, sequence_hashes);
        while (true) {
            const auto range = motifs_counter.get_and_increment_range_info(params.motif_range_size);
            if (range.count() == 0) {
                break;
            }
            motif_finder_gpu_internal(gpu_memory, params, out_motif_weights, range.start, range.count(), thread_id);
        }
    });
}