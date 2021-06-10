#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include <sequence_nums.h>
#include <stat_model.h>

class MarkovStatModel : public StatModel {
public:
    MarkovStatModel(const std::vector<std::string> &sequences,
                    bool complementary,
                    int level,
                    bool use_binom_instead_of_chi2,
                    double binom_correction);

    double motif_probability(uint32_t hash) const override;
    double motif_probability_x4(uint32_t hash) const override;

private:
    void calc_probabilities(uint32_t kmer_length);
    void precalc_probabilities(int symbols_in_key, std::vector<double> &result_probabilites, bool prefix);

private:
    int _level;
    int _kmer_length;
    std::vector<double> _probabilities;
    int _symbols_in_key1;
    int _symbols_in_key2;

    int _key_size1;
    int _key_size2;
    std::vector<double>
        _precalc_prefix_probabilities; // предрассчитанные вероятности для нескольких символов (для префикса)
    std::vector<double>
        _precalc_suffix_probabilities; // предрассчитанные вероятности для нескольких символов (для суффикса)
};
