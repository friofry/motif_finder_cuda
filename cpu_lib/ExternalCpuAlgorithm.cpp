#include "ExternalCpuAlgorithm.h"

#include "CpuCommonAlgorithm.h"
#include <algorithm>

ExternalCpuAlgorithm::ExternalCpuAlgorithm(const std::vector<uint32_t> &motifHashes, const SequenceHashes &sequenceHashes, uint32_t rangeSize)
    : _motifHashes(motifHashes)
    , _sequenceHashes(sequenceHashes)
    , _rangeSize(rangeSize)
{
}

void ExternalCpuAlgorithm::getOccurrenceInRange(uint32_t begin, uint32_t end, std::vector<uint16_t> &outMotifWeights) const
{
    for (uint32_t i = begin; i < end; i++)  {
        outMotifWeights[i] = getMotifOccurrence(_motifHashes[i], _sequenceHashes);
    }
}

uint32_t ExternalCpuAlgorithm::getDesiredRangeSize() const
{
    return _rangeSize;
}

void ExternalCpuAlgorithm::getAllOccurrences(std::vector<uint16_t> &outMotifWeights) const
{
    outMotifWeights.resize(_motifHashes.size(), 0);

    for (uint64_t i = 0; i < TOTAL_MOT; i += _rangeSize)  {
        getOccurrenceInRange(i, std::min<uint32_t >(_motifHashes.size(), i + _rangeSize), outMotifWeights);
    }
}


