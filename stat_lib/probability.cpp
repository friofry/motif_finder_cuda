#include "probability.h"

#include <vector>
#include <algorithm>
#include <cmath>

#include <hash_conversions.h>
#include <letter_conversions.h>

using namespace std;

vector<double> calc_probability_x4(const vector<double> &frequencies_ratio)
{
    uint32_t total_hashes_4 = 1 << 16;
    vector<double> result(total_hashes_4, 0);

    for (uint32_t idx = 0; idx < total_hashes_4; idx++) {
        double prefix_probability = 1.0;
        uint32_t prefix_hash = idx_to_motif_prefix_hash(idx, 4) >> 16;

        for (int pos = 0; pos < 4; pos++) {
            uint32_t cur_symbol = prefix_hash & 0x0000000F;
            if (cur_symbol == 0) {
                prefix_probability = 0.0;
                break;
            }

            double symbol_prob = 0.0;
            for (int nucl = 0; nucl < 4; nucl++) {
                uint8_t nucl_hash = 1 << nucl;

                if (nucl_hash & cur_symbol) {
                    symbol_prob += frequencies_ratio[nucl];
                }
            }

            prefix_probability *= symbol_prob;
            prefix_hash >>= 4;
        }

        result[idx] = prefix_probability;
    }
    return result;
}


vector<uint32_t> calc_frequencies(const std::vector<std::string> &sequences, bool complementary)
{
    vector<uint32_t> result(4, 0);
    for (uint32_t i = 0; i < sequences.size(); ++i) {
        for (char x: sequences[i]) {
            uint8_t nucl_num = to_num(x);
            result[nucl_num]++;
        }
    }

    if (complementary) {
        vector<uint32_t> compl_result(4, 0);
        compl_result[0] = compl_result[0] + compl_result[1];
        compl_result[1] = compl_result[0] + compl_result[1];
        compl_result[2] = compl_result[2] + compl_result[3];
        compl_result[3] = compl_result[2] + compl_result[3];
        return compl_result;
    }
    return result;
}

vector<double> calc_frequencies_ratio(const std::vector<std::string> &sequences, bool complementary)
{
    auto counters = calc_frequencies(sequences, false);

    vector<double> result(4, 0.0);
    uint32_t total = counters[0] + counters[1] + counters[2] + counters[3];
    for (int i = 0; i < 4; i++) {
        result[i] = double(counters[i]) / total;
    }

    if (complementary) {
        vector<double> compl_result(4, 0);
        compl_result[0] = double(counters[0] + counters[1]) / (2.0 * total);
        compl_result[1] = double(counters[0] + counters[1]) / (2.0 * total);
        compl_result[2] = double(counters[2] + counters[3]) / (2.0 * total);
        compl_result[3] = double(counters[2] + counters[3]) / (2.0 * total);
        return compl_result;
    }

    return result;
}


double calc_prob_upper_bound_simple(double external_percentage, unsigned int positions)
{
    return -log(1.0 - external_percentage) / double(positions);
}

double prob_per_sequence_poisson(double prob, unsigned int positions)
{
    // 1 - e ^ (lambda)
    return 1.0 - exp(-prob * positions);
}

double motif_probability(unsigned int motif_hash, const vector<double> &frequencies_ratio)
{
    double res = 1.0;
    for (int i = 0; i < 8; i++) {
        unsigned int cur_symbol = motif_hash & 0xF;
        double symbol_prob = 0.0;
        for (int j = 0; j < 4; j++) {
            if ((1 << j) & cur_symbol) {
                symbol_prob += frequencies_ratio[j];
            }
        }
        res *= symbol_prob;
        motif_hash >>= 4;
    }
    return res;
}


double motif_probability_x4(unsigned int motif_hash, const std::vector<double> &probability_x4)
{
    unsigned int part1 = motif_hash & 0xFFFF;
    unsigned int part2 = (motif_hash >> 16) & 0xFFFF;
    double prob = probability_x4[part1] * probability_x4[part2];
    return prob;
}

/// chi2 criteria
double importance_chi2(
    unsigned int motif_hash,
    unsigned int weight,
    const std::vector<double> &probability_x4,
    unsigned int seq_length,
    unsigned int seq_count)
{
    double position_prob = motif_probability_x4(motif_hash, probability_x4);
    uint32_t expected = prob_per_sequence_poisson(position_prob, seq_length) * seq_count;
    expected = expected ? expected : 1;
    uint32_t x = expected - weight;
    double chi_2 = x * x / expected;
    return chi_2;
}

std::vector<double> calc_probability_x4(
    const std::vector<std::string> &sequences,
    bool search_complementary)
{
    auto frequencies_ratio = calc_frequencies_ratio(sequences, search_complementary);
    return calc_probability_x4(frequencies_ratio);
}
