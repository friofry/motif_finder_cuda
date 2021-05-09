#ifndef MOTIF_FINDER_MOTIFOCCURRENCEEVALUATOR_H
#define MOTIF_FINDER_MOTIFOCCURRENCEEVALUATOR_H

#include <vector>
#include <string>
#include <stdint.h>

#include <IAlgorithm.h>
#include <IAlgorithmFactory.h>
#include <SequenceHashes.h>

// Evaluate motif occurrences in mulltiple threads
// @params outMotifOccurrences - result occurrences
// @params algorithmFactories - creates algorithm instances (for example a hybrid of cpu and gpu algorithms)
// @params sequenceHashes - nucleotide sequences info
// @params motifsCount - number of motifs
void evaluateMotifOccurrencesParallel(std::vector<uint16_t> &outMotifOccurrences,
                                      const std::vector<IAlgorithmFactoryPtr> &algorithmFactories,
                                      const SequenceHashes &sequenceHashes,
                                      uint32_t motifsCount);

#endif //MOTIF_FINDER_MOTIFOCCURRENCEEVALUATOR_H
