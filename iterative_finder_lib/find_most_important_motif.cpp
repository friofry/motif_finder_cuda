#include "find_most_important_motif.h"

#include <thread>
#include <algorithm>

#include <safe_counter.h>
#include <timer.h>

using namespace std;
namespace {

MotifData gather_most_important_motif(const vector<MotifData> &results)
 {
    int max_i = 0;
    for (uint32_t i = 0; i < results.size(); i++) {
        if (results[i].score > results[max_i].score) {
            max_i = i;
        }
    }
    return results[max_i];
}

MotifData find_most_important_motif_in_range(
    uint32_t start_idx,
    uint32_t end_idx,
    const vector<uint16_t> &motif_weights,
    const vector<uint32_t> &motif_hashes,
    const StatModel &stat_model,
    double min_presence,
    bool use_binom_instead_chi2,
    const std::vector<uint16_t> &prev_weights,
    vector<double> &binomial_prob)
{
    double max_score = -100000.0;
    int64_t max_i = -1;

    for (uint64_t i = start_idx; i < end_idx; i++) {
        uint16_t weight = motif_weights[i];
        if (weight == 0) {
            continue;
        }
        uint32_t hash = motif_hashes[i];

        int weight_random = stat_model.get_random_weight(hash);
        if (weight >= min_presence && weight > weight_random) {
            double score = 0.0;
            if (use_binom_instead_chi2) {
                if (prev_weights.size() && weight == prev_weights[i]) {
                    score = binomial_prob[i];
                } else {
                    score = stat_model. binom_by_hash(hash, weight, max_score);
                    binomial_prob[i] = score;
                }
            } else {
                score = stat_model.chi2_by_hash(hash, weight);
            }

            if (score > max_score) {
                max_score = score;
                max_i = i;
            }
        }
    }

    MotifData result;

    if (max_i >= 0) {
        result.hash = motif_hashes[max_i];
        result.weight = motif_weights[max_i];
        result.score = max_score;
    } else {
        result.hash = 0;
        result.weight = 0;
        result.score = -1.0;
    }
    return result;
}

} // namespace

MotifData find_most_important_motif(
    const std::vector<uint16_t> &motif_weights,
    const std::vector<uint32_t> &motif_hashes,
    const StatModel &stat_model,
    double min_presence,
    bool use_binom_instead_chi2,
    const std::vector<uint16_t> &prev_weights,
    vector<double> &binomial_prob)
{
    Timer t("find_most_important_motif");
    t.silence();
    vector<thread> threads;
    SafeCounter counter(motif_weights.size());
    uint64_t total_mots = motif_weights.size();

    uint32_t cpu_thread_count = thread::hardware_concurrency();

    int chunks_per_thread = 2;
    int chunks = chunks_per_thread * cpu_thread_count;
    vector<MotifData> results(chunks);
    int motifs_per_chunk = (total_mots/double(chunks)) + 1;

    for (uint32_t thread_id = 0; thread_id < cpu_thread_count; thread_id++) {
        threads.push_back(thread(
        [&]()
        {
            while (true) {
                auto mot_idx_range_info = counter.get_and_increment_range_info(motifs_per_chunk);

                if (mot_idx_range_info.start >= total_mots) {
                    break;
                }

                auto res = find_most_important_motif_in_range(
                                mot_idx_range_info.start,
                                mot_idx_range_info.end,
                                ref(motif_weights),
                                ref(motif_hashes),
                                ref(stat_model),
                                min_presence,
                                use_binom_instead_chi2,
                                ref(prev_weights),
                                ref(binomial_prob));

                results[mot_idx_range_info.idx] = res;
            }
        }));
    }

    for (uint32_t  i = 0; i < cpu_thread_count; i++) {
        threads[i].join();
    }

    return gather_most_important_motif(results);
}
