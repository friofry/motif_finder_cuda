#include "iterative_finder.h"

#include <algorithm>
#include <fstream>
#include <iostream>

#include <external_generator_cpu.h>
#include <fst_reader.h>
#include <hash_conversions.h>
#include <match_functions.h>
#include <run_parallel.h>
#include <sequence_hashes_utils.h>
#include <sequences_to_hashes.h>

#include "create_stat_model.h"
#include "find_most_important_motif.h"

using namespace std;

ImportantMotifFinder::ImportantMotifFinder(const ArgoCudaParams &params,
                                           const FindOccurrencesAlgorithm &external_algorithm,
                                           const char *output_file)
    : _params(params)
    , _external_algorithm(external_algorithm)
    , _output_file(output_file)
{
    if (params.use_old_motifs_file) {
        throw invalid_argument("use_old_motifs_file is not supported yet.\n");
    }
}

std::vector<uint32_t> ImportantMotifFinder::find()
{
    printf("1. reading sequences\n");
    // 1. Read sequences
    auto sequences = read_fasta(_params.positive_sequences.c_str());
    _sequence_hashes = sequences_to_hashes(sequences, _params.complementary);

    // 2. Init stat model
    printf("2. Init stat model\n");
    _stat_model
        = create_stat_model(_params.use_real_nucl_frequences, _params.markov_level, sequences, _params.complementary);

    // 3. Generate motif hashes
    printf("3. Generate motif hashes\n");
    vector<uint32_t> motif_hashes;
    generate_external_hashes_cpu(_params.max_motif_prob_by_chance, *_stat_model, _params.complementary, motif_hashes);

    // 4. Pre filter hashes
    printf("4. Pre filter hashes\n");
    if (_params.max_motif_score_contrast > 0) {
        try {
            auto contrast_sequences = read_fasta(_params.contrast_sequences.c_str());
            if (!contrast_sequences.empty()) {
                auto contrast_stat_model = create_stat_model(
                    _params.use_real_nucl_frequences, _params.markov_level, contrast_sequences, _params.complementary);
                auto contrast_hashes = sequences_to_hashes(contrast_sequences, _params.complementary);
                vector<uint16_t> contrast_occurrences;
                _external_algorithm(motif_hashes, contrast_hashes, contrast_occurrences);
                filter_hashes_by_contrast(_params.max_motif_score_contrast,
                                          *contrast_stat_model,
                                          contrast_occurrences,
                                          motif_hashes,
                                          _params.use_binom_instead_chi2);
            }
        } catch (const std::invalid_argument &ex) {
            // nop
        }
    }

    // 5. Find and exclude sequences
    printf("5. Find and exclude sequences\n");
    if (_params.chip_seq_percentage > 0.0) {
        find_motifs_iterative(motif_hashes, true);
    }

    // 6. Main search
    printf("6. Main search\n");
    find_motifs_iterative(motif_hashes, false);
    return _found_motifs;
}
void ImportantMotifFinder::find_motifs_iterative(vector<uint32_t> &motif_hashes, bool exclude_sequence)
{
    double min_presence = _params.min_motif_presence * _sequence_hashes.count;
    int minimum_sequences = _params.chip_seq_percentage * _sequence_hashes.count;

    vector<uint16_t> weights;
    if (_found_motifs.size() >= static_cast<uint32_t>(_params.motif_to_find)) {
        return;
    }

    vector<double> binomial_prob;
    vector<uint16_t> prev_weights;
    while (true) {
        if (!_found_motifs.empty()) {
            // 1. Mask found motif hashes in sequences
            printf("1. Mask found motif hashes in sequences\n");
            auto latest_motif_data = _found_motifs_data[_found_motifs.size() - 1];
            motif_hashes[latest_motif_data.index] = 0;
            mask_motif_hashes(_sequence_hashes, latest_motif_data.hash, _params.complementary);

            fill(weights.begin(), weights.end(), 0);
        }
        // 2. Run external algorithm
        printf("2. Run external algorithm\n");
        _external_algorithm(motif_hashes, _sequence_hashes, weights);

        if (_found_motifs_data.empty()) {
            // 3. Exclude motif hashes by chi2
            printf("3. Exclude motif hashes by chi2\n");
            exclude_motifs_by_chi2(motif_hashes, weights);
            if (_params.use_binom_instead_chi2) {
                binomial_prob.resize(motif_hashes.size());
            }
        }

        // 4. Find most significant motif
        printf("4. Find most significant motif\n");
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

        _found_motifs.push_back(max_motif_data.hash);
        _found_motifs_data.push_back(max_motif_data);
        write_results_old();
        if (_found_motifs.size() >= static_cast<uint32_t>(_params.motif_to_find)) {
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

void ImportantMotifFinder::exclude_motifs_by_chi2(std::vector<uint32_t> &motif_hashes,
                                                  std::vector<uint16_t> &weights) const
{
    if (_params.min_motif_chi2 <= 0.0) {
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
            auto chi2 = _stat_model->chi2_by_hash(hash, weight);
            if (chi2 >= _params.min_motif_chi2) {
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

void ImportantMotifFinder::write_results_old()
{
    ofstream f(_output_file.c_str());
    for (uint32_t i = 0; i < _found_motifs_data.size(); i++) {
        const auto &d = _found_motifs_data[i];
        auto rand_w = _stat_model->get_random_weight(d.hash);

        f << hash_to_string(d.hash) << "\t";
        f << int(100 * d.weight / _sequence_hashes.count) << "\t";
        f << int(100 * rand_w / _sequence_hashes.count) << "\t";
        f << int(d.score) << endl;
    }
}

std::vector<uint32_t> find_important_motifs(const ArgoCudaParams &params,
                                            const FindOccurrencesAlgorithm &external_algorithm,
                                            const char *output_file)
{
    ImportantMotifFinder finder(params, external_algorithm, output_file);
    return finder.find();
}