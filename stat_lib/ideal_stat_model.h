#pragma once
#include <vector>

#include <stat_model.h>

// Модель, где все нуклеотиды равновероятны
class IdealStatModel : public StatModel {
public:
    IdealStatModel(const std::vector<std::string> &sequences, bool complementary);

    double motif_probability(uint32_t hash) const override;
    double motif_probability_x4(uint32_t hash) const override;

private:
    // чатсоты для каждого нуклеотида (0.25)
    std::vector<double> _frequencies_ratio;
    // Предрассчитанные вероятности для 4 символов
    std::vector<double> _probability_x4;
};
