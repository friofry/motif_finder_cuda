#ifndef MOTIF_FINDER_CPURANGEFINDER_H
#define MOTIF_FINDER_CPURANGEFINDER_H

#include "i_range_motif_finder.h"
#include "sequence_hashes_params.h"

#include <vector>

class CpuRangeMotifFinder : public IRangeMotifFinder {
public:
    CpuRangeMotifFinder(const SequenceHashesParams &params);
    void get_occurrence_in_range(const std::vector<uint32_t> &motif_hashes, std::vector<uint16_t> &out_motif_weights) const override;

private:
    uint32_t get_occurrence(uint32_t motif_hash) const;

private:
    const SequenceHashesParams &_params;
};

#endif //MOTIF_FINDER_CPURANGEFINDER_H
