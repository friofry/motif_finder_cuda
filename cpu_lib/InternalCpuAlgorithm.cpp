#include "InternalCpuAlgorithm.h"
#include <config.h>
#include <hash_conversions.h>
#include "CpuCommonAlgorithm.h"

InternalCpuAlgorithm::InternalCpuAlgorithm(const SequenceHashes &sequenceHashes, uint32_t rangeSize)
    : _sequenceHashes(sequenceHashes)
    , _rangeSize(rangeSize)
{
}

void InternalCpuAlgorithm::getOccurrenceInRange(uint32_t begin, uint32_t end, std::vector<uint16_t> &outMotifWeights) const
{
    for (uint32_t i = begin; i < end; i++)  {
        uint32_t motifHash = idx_to_hash(i);
        outMotifWeights[i] = getMotifOccurrence(motifHash, _sequenceHashes);
    }
}

uint32_t InternalCpuAlgorithm::getDesiredRangeSize() const
{
    return _rangeSize;
}

void InternalCpuAlgorithm::getAllOccurrences(std::vector<uint16_t> &outMotifWeights) const
{
    outMotifWeights.resize(TOTAL_MOT, 0);
    for (uint64_t i = 0; i < TOTAL_MOT; i += _rangeSize)  {
        getOccurrenceInRange(i, std::min<uint64_t>(TOTAL_MOT, i + _rangeSize), outMotifWeights);
    }
}
