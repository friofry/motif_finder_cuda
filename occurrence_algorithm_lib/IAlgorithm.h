#ifndef MOTIF_FINDER_I_ALGORITHM_H
#define MOTIF_FINDER_I_ALGORITHM_H

#include <stdint.h>
#include <vector>

// Algorithm for evaluating the occurrence of motifs
class IAlgorithm {
public:
    virtual ~IAlgorithm() {}

    // Get motif occurrences and store to outMotifWeights
    // @params begin - first motif index
    // @params end - last motif index + 1 (points after the last motif)
    // @params outMotifWeights - array for results (of at least (end - begin) size)
    virtual void getOccurrenceInRange(uint32_t begin, uint32_t end, std::vector<uint16_t> &outMotifWeights) const = 0;

    // Get all motif occurrences and store to outMotifWeights
    // @params outMotifWeights - array for results (of at least (end - begin) size)
    virtual void getAllOccurrences(std::vector<uint16_t> &outMotifWeights) const = 0;

    // Get desired motifs range
    virtual uint32_t getDesiredRangeSize() const = 0;
};
using IAlgorithmPtr = std::shared_ptr<IAlgorithm>;
#endif // MOTIF_FINDER_I_ALGORITHM_H
