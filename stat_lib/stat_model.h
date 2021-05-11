#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <string>

#include <sequence_nums.h>

// Base class for stat models (abstract)
class StatModel {
public:
    StatModel(const std::vector<std::string> &sequences, bool complementary);

    virtual ~StatModel();

    // Вероятность встретить мотив по случайным причинам
    virtual double motif_probability(uint32_t hash) const = 0;

    // оптимизированная версия
    virtual double motif_probability_x4(uint32_t hash) const = 0;

    // Хи квадрат (без эвристик)
    double chi2(double weight, double random_weight) const;

    // Возвращает ожидаемую встречаемость (без эвристик)
    double get_random_weight(uint32_t hash) const;

    double chi2_by_hash(uint32_t hash, uint16_t weight) const;

    int get_avg_hashes_per_sequence() const;

    // Биномиальный критерий (Бернулли)
    double binom_by_hash(uint32_t hash, uint16_t weight, double max_score = 0.0) const;

protected:
    std::vector<std::string> _sequences;
    int _hashes_per_sequence { 0 };
    SequenceNums _sequence_nums;
    bool _complementary;
    std::vector<double> _log_sums;
};

typedef std::shared_ptr<StatModel> StatModelPtr;