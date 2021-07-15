#include "external_generator_cpu.h"

#include <cmath>
#include <iostream>
#include <thread>

#include <config.h>
#include <hash_conversions.h>
#include <hash_conversions_x4.h>
#include <timer.h>

#include "probability.h"

using namespace std;

void generate_external_indexes_cpu(double max_motif_prob_simple,
                                   const vector<double> &probability_x4,
                                   std::vector<uint32_t> &result)
{
    result.reserve(TOTAL_MOT);

    // precalculate hashes for 4-digit indexes
    auto hash_x4 = calc_hash_x4();

    for (uint32_t i = 0; i < TOTAL_MOT; i++) {
        uint32_t hash = idx_to_hash_x4(i, hash_x4);
        bool important = motif_probability_x4(hash, probability_x4) <= max_motif_prob_simple;
        if (important) {
            result.push_back(i);
        }
    }

    result.shrink_to_fit();
}

void generate_external_indexes_cpu(double external_percentage,
                                   const std::vector<std::string> &sequences,
                                   bool search_complementary,
                                   std::vector<uint32_t> &result)
{
    if (sequences.empty()) {
        cout << "Error, empty sequences list" << endl;
        return;
    }

    auto probability_x4 = calc_probability_x4(sequences, search_complementary);

    uint32_t positions = sequences[0].length() - MOTIV_LEN + 1;
    double max_motif_prob_simple = calc_prob_upper_bound_simple(external_percentage, positions);

    return generate_external_indexes_cpu(max_motif_prob_simple, probability_x4, result);
}

namespace {

double calculate_border_motif_probability(double max_motif_prob_by_chance, int number_of_hashes)
{
    return 1.0 - exp(log(1.0 - max_motif_prob_by_chance) / number_of_hashes);
}

void generate_hashes_thread(vector<uint32_t> &result_hashes,
                            const StatModel &stat_model,
                            bool complementary,
                            double probability_border,
                            uint32_t start,
                            uint32_t end)
{
    auto hash_x4 = calc_hash_x4();
    for (uint32_t i = start; i < end; i++) {
        uint32_t hash = idx_to_hash_x4(i, hash_x4);

        double probability = stat_model.motif_probability_x4(hash);
        if (probability > probability_border) {
            continue;
        }
        if (complementary) {
            uint32_t compl_hash = to_compl_hash_reverse(hash);
            // the complementary motif has been already processed
            if (compl_hash < hash) {
                continue;
            }
        }

        result_hashes.push_back(hash);
    }
    result_hashes.shrink_to_fit();
}

void generate_hashes(vector<uint32_t> &result_hashes,
                     const StatModel &stat_model,
                     bool complementary,
                     double probability_border)
{
    vector<thread> threads;
    int threads_count = thread::hardware_concurrency();
    vector<vector<uint32_t>> thread_results(threads_count);
    for (auto &t : thread_results) {
        t.reserve(TOTAL_MOT / threads_count + 1);
    }

    uint32_t chunk_size = TOTAL_MOT / (threads_count - 1);

    for (int i = 0; i < threads_count; i++) {
        uint32_t start = i * chunk_size;
        uint32_t end = start + chunk_size;
        if (end > TOTAL_MOT) {
            end = TOTAL_MOT;
        }

        threads.push_back(std::thread(generate_hashes_thread,
                                      ref(thread_results[i]),
                                      ref(stat_model),
                                      complementary,
                                      probability_border,
                                      start,
                                      end));
    }

    for (int i = 0; i < threads_count; i++) {
        threads[i].join();
    }

    result_hashes.clear();
    uint32_t sz = 0;
    for (const auto &res : thread_results) {
        sz += res.size();
    }
    result_hashes.reserve(sz);

    for (auto &res : thread_results) {
        result_hashes.insert(result_hashes.end(), res.begin(), res.end());
        res.clear();
    }
    result_hashes.shrink_to_fit();
}

void filter_hashes_thread(vector<uint32_t> &result_hashes,
                          const StatModel &stat_model,
                          const vector<uint16_t> &weights,
                          const vector<uint32_t> &motifs_hashes,
                          double max_score,
                          uint32_t start,
                          uint32_t end)
{
    for (uint32_t i = start; i < end; i++) {
        uint16_t weight = weights[i];
        uint32_t hash = motifs_hashes[i];
        double score = stat_model.score(hash, weight);

        if (score < max_score) {
            result_hashes.push_back(hash);
        }
    }
}

void filter_hashes(double max_score,
                   const StatModel &stat_model,
                   const vector<uint16_t> &weights,
                   vector<uint32_t> &motifs_hashes)
{
    vector<thread> threads;
    int threads_count = thread::hardware_concurrency();
    vector<vector<uint32_t>> thread_results(threads_count);
    for (auto &t : thread_results) {
        t.reserve(MOT_PER_CHUNK);
    }

    uint32_t chunk_size = motifs_hashes.size() / (threads_count - 1);

    for (int i = 0; i < threads_count; i++) {
        uint32_t start = i * chunk_size;
        uint32_t end = start + chunk_size;
        if (end > motifs_hashes.size()) {
            end = motifs_hashes.size();
        }

        threads.push_back(std::thread(filter_hashes_thread,
                                      ref(thread_results[i]),
                                      ref(stat_model),
                                      ref(weights),
                                      ref(motifs_hashes),
                                      max_score,
                                      start,
                                      end));
    }

    for (int i = 0; i < threads_count; i++) {
        threads[i].join();
    }

    motifs_hashes.clear();
    uint32_t sz = 0;
    for (const auto &res : thread_results) {
        sz += res.size();
    }
    motifs_hashes.reserve(sz);

    for (const auto &res : thread_results) {
        motifs_hashes.insert(motifs_hashes.end(), res.begin(), res.end());
    }
    motifs_hashes.shrink_to_fit();
}

} // namespace

void generate_external_hashes_cpu(double max_motif_prob_by_chance,
                                  const StatModel &stat_model,
                                  bool complementary,
                                  std::vector<uint32_t> &result_hashes)
{
    Timer t("generate_external_hashes_cpu");
    t.silence();
    int hashes_per_sequence = stat_model.get_avg_hashes_per_sequence();

    // Boundary value of the probability of occurrence of the motif in the position
    double probability_border = calculate_border_motif_probability(max_motif_prob_by_chance, hashes_per_sequence);
    // printf("prob border %f %u: %f\n", max_motif_prob_by_chance, hashes_per_sequence, probability_border);

    generate_hashes(result_hashes, stat_model, complementary, probability_border);
}

void filter_hashes_by_contrast(double max_score,
                               const StatModel &stat_model,
                               const vector<uint16_t> &weights,
                               vector<uint32_t> &motifs_hashes)
{
    Timer t("filter_indexes_by_contrast");
    t.silence();
    filter_hashes(max_score, stat_model, weights, motifs_hashes);
}

void generate_all_hashes(vector<uint32_t> &motif_hashes, uint32_t count)
{
    vector<thread> threads;
    int threads_count = thread::hardware_concurrency();

    motif_hashes.resize(count);
    uint32_t chunk_size = count / (threads_count - 1);
    auto hash_x4 = calc_hash_x4();

    for (int i = 0; i < threads_count; i++) {
        uint32_t start = i * chunk_size;
        uint32_t end = start + chunk_size;
        if (end > count) {
            end = count;
        }
        threads.push_back(std::thread([start, end, &motif_hashes, &hash_x4]() {
            for (uint32_t i = start; i < end; i++) {
                motif_hashes[i] = idx_to_hash_x4(i, hash_x4);
            }
        }));
    }

    for (int i = 0; i < threads_count; i++) {
        threads[i].join();
    }
}
