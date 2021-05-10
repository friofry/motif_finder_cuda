#include "cpu_common_algorithm.h"

uint32_t get_motif_occurrence(uint32_t motif_hash, const SequenceHashes &params)
{
    uint32_t result = 0;

    for (uint32_t i = 0; i < params.count; i++) {
        for (uint32_t pos = 0; pos < params.lengths[i]; pos++) {
            uint32_t seq_hash = params.hashes[params.seq_begins[i] + pos];
            if ((seq_hash & motif_hash) == seq_hash) {
                result++;
                break;
            }
        }
    }
    return result;
}