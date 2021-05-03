#include "cpu_range_motif_finder.h"

CpuRangeMotifFinder::CpuRangeMotifFinder(const SequenceHashesParams &params)
    :_params(params)
{
}

void CpuRangeMotifFinder::get_occurrence_in_range(const std::vector<uint32_t> &motif_hashes, std::vector<uint16_t> &out_motif_weights) const
{
    for (uint32_t i = 0; i < motif_hashes.size(); i++)  {
        out_motif_weights[i] = get_occurrence(motif_hashes[i]);
    }
}

uint32_t CpuRangeMotifFinder::get_occurrence(uint32_t motif_hash) const
{
    uint32_t result = 0;

    for (uint32_t i = 0; i < _params.count; i++) {
        uint32_t seq_len = _params.lengths[i];
        uint32_t seq_start = _params.seq_begins[i];

        for (uint32_t pos = 0; pos < seq_len; pos++) {
            uint32_t seq_hash = _params.hashes[seq_start + pos];
            if ((seq_hash & motif_hash) == seq_hash) {
                result++;
                break;
            }
        }
    }
    return result;
}