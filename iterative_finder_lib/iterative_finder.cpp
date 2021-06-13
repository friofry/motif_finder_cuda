#include "iterative_finder.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>

#include <external_generator_cpu.h>
#include <fst_reader.h>
#include <hash_conversions.h>
#include <match_functions.h>
#include <run_parallel.h>
#include <sequence_hashes_utils.h>
#include <sequences_to_hashes.h>

#include "create_stat_model.h"
#include "filter_shifted_motifs.h"
#include "find_most_important_motif.h"
#include "motif_data_utils.h"

using namespace std;

ImportantMotifFinder::ImportantMotifFinder(const ArgoCudaParams &params,
                                           const FindOccurrencesAlgorithm &external_algorithm)
    : _params(params)
    , _external_algorithm(external_algorithm)
    , _output_file(params.output_file)
{
    if (_output_file.empty()) {
        _output_file = "a.txt";
    }

    if (params.use_old_motifs_file) {
        throw invalid_argument("use_old_motifs_file is not supported yet.\n");
    }

    // 1. Read sequences
    printf("Init: Read sequences\n");
    auto sequences = read_fasta(_params.positive_sequences.c_str());
    _sequence_hashes = sequences_to_hashes(sequences, _params.complementary);

    // 2. Init stat model
    double correction = _params.bonferroni_correction ? -log10(TOTAL_MOT) : 0.0;
    printf("Init: Create stat model. Real freqs: %d, Markov level: %d, Sequence count: %lu, Complementary: %d, Use binom score: %d, Score correction: %f\n",
           _params.use_real_nucl_frequences,
           _params.markov_level,
           sequences.size(),
           _params.complementary,
           _params.use_binom_instead_chi2,
           correction);

    _stat_model = create_stat_model(_params.use_real_nucl_frequences,
                                    _params.markov_level,
                                    sequences,
                                    _params.complementary,
                                    _params.use_binom_instead_chi2,
                                    correction);
}

std::vector<uint32_t> ImportantMotifFinder::find()
{
    // 3. Generate motif hashes
    printf("Init: Generating motif hashes...\n");
    vector<uint32_t> motif_hashes;
    generate_external_hashes_cpu(_params.max_motif_prob_by_chance, *_stat_model, _params.complementary, motif_hashes);
    printf("Init: Motif hashes %lu\n", motif_hashes.size());
    // 4. Pre filter hashes
    if (_params.max_motif_score_contrast > 0) {
        try {
            auto contrast_sequences = read_fasta(_params.contrast_sequences.c_str());
            if (!contrast_sequences.empty()) {
                double correction = _params.bonferroni_correction ? -log10(TOTAL_MOT) : 0.0;
                auto contrast_stat_model = create_stat_model(_params.use_real_nucl_frequences,
                                                             _params.markov_level,
                                                             contrast_sequences,
                                                             _params.complementary,
                                                             _params.use_binom_instead_chi2,
                                                             correction);
                auto contrast_hashes = sequences_to_hashes(contrast_sequences, _params.complementary);
                vector<uint16_t> contrast_occurrences;
                printf("Init: Finding occurrences in the contrast sequence set...\n");
                _external_algorithm(motif_hashes, contrast_hashes, contrast_occurrences);

                auto hashes_size = motif_hashes.size();
                filter_hashes_by_contrast(
                    _params.max_motif_score_contrast, *contrast_stat_model, contrast_occurrences, motif_hashes);
                printf("Init: Filtered hashes %lu\n", hashes_size - motif_hashes.size());
            }
        } catch (const std::invalid_argument &ex) {
            // nop
        }
    }

    // 5. Find and exclude sequences
    printf("Init: Find and exclude sequences\n");
    if (_params.chip_seq_percentage > 0.0) {
        find_motifs_iterative(motif_hashes, true);
    }

    // 6. Main search
    printf("Main search:\n\n");
    find_motifs_iterative(motif_hashes, false);
    return get_hashes(_found_motifs_data);
}
void ImportantMotifFinder::find_motifs_iterative(vector<uint32_t> &motif_hashes, bool exclude_sequence)
{
    double min_presence = _params.min_motif_presence * _sequence_hashes.count;
    int minimum_sequences = _params.chip_seq_percentage * _sequence_hashes.count;

    vector<uint16_t> weights;
    if (_found_motifs_data.size() >= static_cast<uint32_t>(_params.motif_to_find)) {
        return;
    }

    vector<double> binomial_prob;
    vector<uint16_t> prev_weights;
    while (true) {
        if (!_found_motifs_data.empty()) {
            // 1. Mask found motif hashes in sequences
            auto latest_motif_data = _found_motifs_data[_found_motifs_data.size() - 1];
            printf("1. Mask found motif hashes in sequences %s\n", hash_to_string(latest_motif_data.hash).c_str());
            motif_hashes[latest_motif_data.index] = 0;
            remove_motif_hashes(_sequence_hashes, latest_motif_data.hash, _params.complementary);
            auto s = hashes_to_sequences(_sequence_hashes, _params.complementary);
            fill(weights.begin(), weights.end(), 0);
        }
        // 2. Run external algorithm
        printf("2. Running external algorithm %lu...\n", motif_hashes.size());
        _external_algorithm(motif_hashes, _sequence_hashes, weights);

        if (_found_motifs_data.empty()) {
            // 3. Exclude motif hashes by score
            printf("3. Excluding motif hashes by score...\n");
            exclude_motifs_by_score(motif_hashes, weights);
            printf("3. Motifs left: %lu\n", motif_hashes.size());
            if (_params.use_binom_instead_chi2) {
                binomial_prob.resize(motif_hashes.size());
            }
        }

        // 4. Find most significant motif
        printf("4. Find motif with max score\n");
        auto max_motif_data = find_most_important_motif(weights,
                                                        motif_hashes,
                                                        *_stat_model,
                                                        min_presence,
                                                        _params.use_binom_instead_chi2,
                                                        prev_weights,
                                                        binomial_prob);
        if (_params.use_binom_instead_chi2) {
            prev_weights.resize(weights.size());
            std::copy(weights.begin(), weights.end(), prev_weights.begin());
        }

        if (max_motif_data.score < _params.min_motif_score) {
            break;
        }

        _found_motifs_data.push_back(max_motif_data);
        write_results();
        if (_found_motifs_data.size() >= static_cast<uint32_t>(_params.motif_to_find)) {
            break;
        }

        if (exclude_sequence) {
            int sequences_left = nonempty_sequence_count(_sequence_hashes);
            if (sequences_left < minimum_sequences || sequences_left == 0) {
                break;
            }
        }
    }
}

void ImportantMotifFinder::exclude_motifs_by_score(std::vector<uint32_t> &motif_hashes,
                                                   std::vector<uint16_t> &weights) const
{
    if (_params.min_motif_score <= 0.0) {
        return;
    }
    double min_presence = _params.min_motif_presence * _sequence_hashes.count;
    std::vector<uint32_t> result_ids;
    result_ids.reserve(motif_hashes.size());

    for (uint32_t i = 0; i < motif_hashes.size(); i++) {
        auto weight = weights[i];
        auto hash = motif_hashes[i];
        int weight_random = _stat_model->get_random_weight(hash);
        if (weight >= min_presence && weight > weight_random) {
            auto score = _stat_model->score(hash, weight);
            if (score >= _params.min_motif_score) {
                result_ids.push_back(i);
            }
        }
    }
    printf("Ratio: %f \n", (100.0 * result_ids.size() / motif_hashes.size()));
    std::vector<uint32_t> motif_hashes_new(result_ids.size());
    for (uint32_t i = 0; i < result_ids.size(); i++) {
        auto idx = result_ids[i];
        motif_hashes_new[i] = motif_hashes[idx];
    }
    motif_hashes.swap(motif_hashes_new);

    std::vector<uint16_t> weights_new(result_ids.size());
    for (uint32_t i = 0; i < result_ids.size(); i++) {
        auto idx = result_ids[i];
        weights_new[i] = weights[idx];
    }
    weights.swap(weights_new);
}

void ImportantMotifFinder::write_results()
{
    ofstream f(_output_file.c_str());
    auto hashes = get_hashes(_found_motifs_data);
    for (uint32_t i = 0; i < _found_motifs_data.size(); i++) {
        if (_params.skip_shifted_results && hash_matches(_found_motifs_data[i].hash, hashes, _params.complementary)) {
            continue;
        }
        const auto &d = _found_motifs_data[i];
        auto rand_w = _stat_model->get_random_weight(d.hash);
        f << hash_to_string(d.hash) << "\t";
        f << int(100 * d.weight / _sequence_hashes.count) << "\t";
        f << int(100 * rand_w / _sequence_hashes.count) << "\t5";
        f << (_params.int_results ? int(d.score) : d.score) << endl;
    }
}

std::vector<uint32_t> find_important_motifs(const ArgoCudaParams &params,
                                            const FindOccurrencesAlgorithm &external_algorithm)
{
    ImportantMotifFinder finder(params, external_algorithm);
    return finder.find();
}