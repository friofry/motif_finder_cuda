#ifndef MOTIF_FINDER_IALGORITHM_H
#define MOTIF_FINDER_IALGORITHM_H

#include <memory>
#include <IAlgorithm.h>
#include <SequenceHashesParams.h>

class IAlgorithmFactory {
public:
    virtual ~IAlgorithmFactory() {}

    // Get occurrence evaluator for specific thread
    virtual IAlgorithmPtr create(int thread_idx, const SequenceHashesParams &hashesInfo) const = 0;

    // Get number of motifs to process
    virtual uint32_t motifsCount() const = 0;

    // Get number of threads
    virtual uint32_t threadsCount() const = 0;
};

using IAlgorithmPtr = std::shared_ptr<IAlgorithm>;
#endif //MOTIF_FINDER_IALGORITHM_H
