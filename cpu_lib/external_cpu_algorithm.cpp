#include "external_cpu_algorithm.h"

#include <algorithm>
#include <thread>

#include <config.h>
#include <run_parallel.h>
#include <safe_counter.h>

#include "cpu_common_algorithm.h"

void external_cpu_algorithm(const std::vector<uint32_t> &motif_hashes,
                            const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights,
                            uint32_t threads)
{
    out_motif_weights.resize(motif_hashes.size(), 0);

    threads = (threads > 0) ? threads : std::thread::hardware_concurrency();

    SafeCounter motifs_counter(motif_hashes.size());

    run_parallel(threads, [&](uint32_t ) {
        while (true) {
            auto range = motifs_counter.get_and_increment_range_info(MOT_PER_CHUNK / 100);
            if (range.count() == 0) {
                break;
            }
            for (uint32_t i = range.start; i < range.end; i++) {
                out_motif_weights[i] = get_motif_occurrence(motif_hashes[i], sequence_hashes);
            }
        }
    });
}
