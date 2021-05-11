#pragma once
#include <vector>
#include <string>

#include <stat_model.h>

class BernulliStatModel : public StatModel {
public:
    BernulliStatModel(const std::vector<std::string> &sequences, bool complementary);

    double motif_probability(uint32_t hash) const override;
    double motif_probability_x4(uint32_t hash) const override;

private:
    std::vector<double> _frequencies_ratio;
    // Предрассчитанные вероятности для 4 символов
    std::vector<double> _probability_x4;
};
