#include "stat_model.h"
#include <algorithm>
#include <cmath>

#include <config.h>
#include <fst_reader.h>
#include <sequences_to_nums.h>

using namespace std;

StatModel::StatModel(const vector<string> &sequences, bool complementary, bool use_binom_instead_of_chi2, double binom_correction)
    : _sequences(sequences)
    , _complementary(complementary)
    , _use_binom_instead_of_chi2(use_binom_instead_of_chi2)
    , _binom_correction(binom_correction)
{
    _sequence_nums = sequences_to_nums(sequences);

    uint32_t average_sequence_length = static_cast<uint32_t>(_sequence_nums.nums.size()) / _sequence_nums.count;
    _hashes_per_sequence = (average_sequence_length - MOTIV_LEN + 1) * (_complementary ? 2 : 1);

    int n = _sequence_nums.count;

    _log_sums.resize(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        _log_sums[i] = _log_sums[i - 1] + log10(double(n - i + 1) / double(i));
    }
}

StatModel::~StatModel() {}

int StatModel::get_avg_hashes_per_sequence() const
{
    return _hashes_per_sequence;
}

double StatModel::get_random_weight(uint32_t hash) const
{
    double p = motif_probability_x4(hash);
    // Expected occurrence of motif in the sequences (count)
    return _sequence_nums.count * (1.0 - exp(_hashes_per_sequence * log(1.0 - p)));
}

double StatModel::chi2(double random_weight, double weight) const
{
    double chi2 = weight - random_weight;
    chi2 = chi2 * chi2 / random_weight;
    return chi2;
}

double StatModel::chi2_by_hash(uint32_t hash, uint16_t weight) const
{
    double random_weight = get_random_weight(hash);
    return chi2(random_weight, weight);
}

double StatModel::binom_by_hash(uint32_t hash, uint16_t k, double max_score) const
{
    double p_in_pos = motif_probability_x4(hash);
    double p = 1.0 - exp(_hashes_per_sequence * log(1.0 - p_in_pos));
    double q = 1.0 - p;
    int n = _sequence_nums.count;
    double lpq = log10(p / q);
    if (k == 0) {
        return 0;
    }

    double a = n * -log10(q);
    double base_log = a - lpq * k - _log_sums[k];

    // 1. find base
    for (int i = k; i <= n; i++) {
        double cur_log = a - lpq * i - _log_sums[i];
        if (cur_log < base_log) {
            base_log = cur_log;
        }
    }

    if (base_log + StatModel::score_difference_drop < max_score) {
        return base_log;
    }
    // 2. find result
    int count = 0;
    double digits = 5;
    double res = 0;
    int max_components = 1000;
    for (int i = k; i <= n; i++) {
        double cur_log = a - lpq * i - _log_sums[i];
        double diff = cur_log - base_log;
        if (diff <= digits) {
            res += pow(0.1, diff);
            count++;
            if (count >= max_components) {
                break;
            }
        }
    }
    return -log10(res) + base_log + _binom_correction;
}

double StatModel::score(uint32_t hash, uint16_t weight, double max_score) const
{
    if (_use_binom_instead_of_chi2) {
        return binom_by_hash(hash, weight, max_score);
    } else {
        return chi2_by_hash(hash, weight);
    }
}
