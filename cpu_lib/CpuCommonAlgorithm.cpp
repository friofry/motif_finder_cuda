#include "CpuCommonAlgorithm.h"

uint32_t getMotifOccurrence(uint32_t motifHash, const SequenceHashes &params)
{
    uint32_t result = 0;

    for (uint32_t i = 0; i < params.count; i++) {
        for (uint32_t pos = 0; pos < params.lengths[i]; pos++) {
            uint32_t seqHash = params.hashes[params.seq_begins[i] + pos];
            if ((seqHash & motifHash) == seqHash) {
                result++;
                break;
            }
        }
    }
    return result;
}