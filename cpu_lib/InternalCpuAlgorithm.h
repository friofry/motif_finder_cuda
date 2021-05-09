#ifndef MOTIF_FINDER_INTERNALCPUALGORITHM_H
#define MOTIF_FINDER_INTERNALCPUALGORITHM_H

#include <IAlgorithm.h>
#include <SequenceHashes.h>
#include <config.h>

class InternalCpuAlgorithm : public IAlgorithm {
public:
    InternalCpuAlgorithm(const SequenceHashes &sequenceHashes, uint32_t rangeSize = MOT_PER_CHUNK/20);

    // Get motif occurrences and store to outMotifWeights
    // @params begin - first motif index
    // @params end - last motif index + 1 (points after the last motif)
    // @params outMotifWeights - array for results (of at least (end - begin) size)
    void getOccurrenceInRange(uint32_t begin, uint32_t end, std::vector<uint16_t> &outMotifWeights) const override;

    // Get desired motifs range
    uint32_t getDesiredRangeSize() const override;

    // Get all motif occurrences and store to outMotifWeights
    // @params outMotifWeights - array for results (of at least (end - begin) size)
    void getAllOccurrences(std::vector<uint16_t> &outMotifWeights) const override;

private:
    const SequenceHashes &_sequenceHashes;
    uint32_t _rangeSize;
};

#endif //MOTIF_FINDER_INTERNALCPUALGORITHM_H
