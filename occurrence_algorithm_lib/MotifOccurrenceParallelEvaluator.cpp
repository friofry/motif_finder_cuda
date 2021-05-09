#include "MotifOccurrenceParallelEvaluator.h"

#include <thread>
#include <safe_counter.h>

using namespace std;

namespace {
void evaluateThreadFunc(std::vector<uint16_t> &outMotifOccurrences,
                        const IAlgorithmFactory &factory,
                        const SequenceHashes &hashesInfo,
                        uint32_t threadId,
                        SafeCounter &rangeCounter)
{
    auto evaluator = factory.create(threadId, hashesInfo);

    const uint32_t motifRangeSize = evaluator.getDesiredRangeSize();
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
                                      const IAlgorithmFactory &algorithmFactory,
                                      const SequenceHashesParams &sequenceHashes)
{
    outMotifOccurrences.resize(algorithmFactory.motifsCount(), 0);

    // Run in parallel
    SafeCounter rangeCounter(params.motif_hashes.size());

    vector<thread> threads;
    int thread_count = algorithmFactory.threadsCount() ? algorithmFactory.threadsCount() : 1;
    for (int i = 0; i < thread_count; i++) {
        threads.push_back(std::thread(evaluateThreadFunc,
                                      ref(outMotifOccurrences),
                                      ref(algorithmFactory),
                                      ref(sequenceHashes),
                                      i,
                                      ref(motifs_counter)));
    }

    // Wait for finished
    for (int i = 0; i < thread_count; i++) {
        threads[i].join();
    }
}