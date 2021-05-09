#include "MotifOccurrenceParallelEvaluator.h"

#include <thread>
#include <safe_counter.h>

using namespace std;

namespace {
void evaluateThreadFunc(std::vector<uint16_t> &outMotifOccurrences,
                        const IAlgorithmFactoryPtr &factory,
                        const SequenceHashes &hashesInfo,
                        uint32_t threadId,
                        SafeCounter &rangeCounter)
{
    auto evaluator = factory->create(threadId, hashesInfo);

    const uint32_t motifRangeSize = evaluator->getDesiredRangeSize();
    std::vector<uint16_t> occurrences(motifRangeSize, 0);
    while (true) {
        auto range_info = rangeCounter.get_and_increment_range_info(motifRangeSize);
        if (range_info.start >= range_info.end) {
            break;
        }

        uint32_t count = range_info.end - range_info.start;
        occurrences.resize(count, 0);
        evaluator->getOccurrenceInRange(range_info.start, range_info.end, occurrences);
        std::copy(occurrences.begin(), occurrences.end(), outMotifOccurrences.begin() + range_info.start);
    }
}
} // namespace

void evaluateMotifOccurrencesParallel(std::vector<uint16_t> &outMotifOccurrences,
                                      const std::vector<IAlgorithmFactoryPtr> &algorithmFactories,
                                      const SequenceHashes &sequenceHashes,
                                      uint32_t motifsCount)
{
    outMotifOccurrences.resize(motifsCount, 0);

    // Run in parallel
    SafeCounter rangeCounter(motifsCount);
    vector<thread> threads;

    for (const auto &factory: algorithmFactories) {
        uint32_t threadCount = factory->threadsCount() ? factory->threadsCount() : 1;
        for (uint32_t i = 0; i < threadCount; i++) {
            threads.push_back(std::thread(evaluateThreadFunc,
                                          ref(outMotifOccurrences),
                                          ref(factory),
                                          ref(sequenceHashes),
                                          i,
                                          ref(rangeCounter)));
        }
    }

    // Wait for finished
    for (uint32_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }
}