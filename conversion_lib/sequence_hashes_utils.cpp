#include "sequence_hashes_utils.h"

#include <cstdint>
#include <set>
#include <thread>
#include <vector>

#include <match_functions.h>
#include <run_parallel.h>
#include <safe_counter.h>

using namespace std;

int nonempty_sequence_count(const SequenceHashes &sequences_hashes)
{
    set<uint32_t> ids;
    int count = 0;
    for (uint32_t i = 0; i < sequences_hashes.count; i++) {
        uint32_t start = sequences_hashes.seq_begins[i];
        uint32_t len = sequences_hashes.lengths[i];
        bool all_zeroes = true;
        for (uint32_t j = 0; j < len; j++) {
            if (sequences_hashes.hashes[start + j] != 0) {
                all_zeroes = false;
                break;
            }
        }
        if (!all_zeroes) {
            count++;
        }
    }
    return count;
}

void mask_motif_hashes(SequenceHashes &sequences_hashes, uint32_t latest_motif, bool complementary)
{
    vector<thread> threads;
    SequenceHashes result;

    SafeCounter counter(sequences_hashes.hashes.size());

    int thread_count = thread::hardware_concurrency();
    uint32_t chunk_size = sequences_hashes.hashes.size() / (2 * thread_count);

    auto fn = [&counter, &sequences_hashes, complementary, latest_motif, chunk_size](uint32_t) {
        while (true) {
            auto range = counter.get_and_increment_range(chunk_size);
            if (range.first >= sequences_hashes.hashes.size()) {
                break;
            }

            for (uint32_t i = range.first; i < range.second; i++) {
                auto hash = sequences_hashes.hashes[i];
                if (hash_matches_motif(hash, latest_motif, complementary)) {
                    sequences_hashes.hashes[i] = 0;
                }
            }
        }
    };

    run_parallel(thread_count, fn);
}

void remove_motif_hashes(SequenceHashes &sequences_hashes, uint32_t latest_motif, bool complementary)
{
    vector<thread> threads;
    SequenceHashes result;
    result.hashes.reserve(sequences_hashes.hashes.size());

    for (uint32_t i = 0; i < sequences_hashes.count; i++) {
        uint32_t len = 0;
        result.seq_begins.push_back(result.hashes.size());
        for (uint32_t j = 0; j < sequences_hashes.lengths[i]; j++) {
            uint32_t hash = sequences_hashes.hashes[j + sequences_hashes.seq_begins[i]];
            if (hash_matches_motif(hash, latest_motif, complementary)) {
                continue;
            }
            result.hashes.push_back(hash);
            len++;
        }
        if (len == 0) {
            result.seq_begins.pop_back();
        }
        result.lengths.push_back(len);
    }
    result.count = result.seq_begins.size();
    std::swap(result, sequences_hashes);
}
