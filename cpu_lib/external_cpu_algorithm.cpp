#include "external_cpu_algorithm.h"

#include "cpu_common_algorithm.h"
#include <algorithm>

void external_cpu_algorithm(const std::vector<uint32_t> &motif_hashes,
                            const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights)
{
    out_motif_weights.resize(motif_hashes.size(), 0);

    for (uint64_t i = 0; i < motif_hashes.size(); i++)  {
        out_motif_weights[i] = get_motif_occurrence(motif_hashes[i], sequence_hashes);
    }
}
