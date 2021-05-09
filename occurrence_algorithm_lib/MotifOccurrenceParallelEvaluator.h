#ifndef MOTIF_FINDER_MOTIFOCCURRENCEEVALUATOR_H
#define MOTIF_FINDER_MOTIFOCCURRENCEEVALUATOR_H

#include <vector>
#include <string>
#include <stdint.h>

#include <IAlgorithm.h>
#include <IAlgorithmFactory.h>
#include <SequenceHashesParams.h>

// Evaluate motif occurrences in mulltiple threads
// @params outMotifOccurrences - result occurrences
// @params algorithmFactory - create algorithm instances
// @params sequenceHashes - nucleotide sequences info
void evaluateMotifOccurrencesParallel(std::vector<uint16_t> &outMotifOccurrences,
                                      const IAlgorithmFactory &algorithmFactory,
                                      const SequenceHashesParams &sequenceHashes);

#endif //MOTIF_FINDER_MOTIFOCCURRENCEEVALUATOR_H
