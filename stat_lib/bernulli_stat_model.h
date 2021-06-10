#pragma once
#include <string>
#include <vector>

#include <stat_model.h>

class BernulliStatModel : public StatModel {
public:
    BernulliStatModel(const std::vector<std::string> &sequences, bool complementary, bool use_binom_instead_of_chi2, double binom_correction);

    double motif_probability(uint32_t hash) const override;
    double motif_probability_x4(uint32_t hash) const override;

private:
    std::vector<double> _frequencies_ratio;
    // Предрассчитанные вероятности для 4 символов
    std::vector<double> _probability_x4;
};
