#include "internal_cpu_algorithm.h"
#include <config.h>
#include <hash_conversions.h>
#include "cpu_common_algorithm.h"

void internal_cpu_algorithm(const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights)
{
    out_motif_weights.resize(TOTAL_MOT, 0);
    for (uint64_t i = 0; i < TOTAL_MOT; i++)  {
        out_motif_weights[i] = get_motif_occurrence(idx_to_hash(i), sequence_hashes);
    }
}
