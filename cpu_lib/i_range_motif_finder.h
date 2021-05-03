#ifndef MOTIF_FINDER_I_RANGE_MOTIF_FINDER_H
#define MOTIF_FINDER_I_RANGE_MOTIF_FINDER_H

#include <stdint.h>
#include <vector>

class IRangeMotifFinder {
public:
    virtual ~IRangeMotifFinder() {}

    // Find occurrences for motifs in the given array and write results to output
    // @param motif_hashes - list of motifs
    // @param out_motif_weights - result occurrences
    virtual void get_occurrence_in_range(const std::vector<uint32_t> &motif_hashes, std::vector<uint16_t> &out_motif_weights) const = 0;
};

#endif //MOTIF_FINDER_I_RANGE_MOTIF_FINDER_H
