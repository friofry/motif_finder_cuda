#pragma once
#include <vector>
#include <string>
#include <cstdint>

double prob_per_sequence_poisson(double prob, unsigned int positions);
double motif_probability(uint32_t motif_hash, const std::vector<double> &frequencies_ratio);

// Оптимизированная функция рассчета вероятности встретить по случайным причинам
double motif_probability_x4(uint32_t motif_hash, const std::vector<double> &probability_x4);

double importance_chi2(
    uint32_t motif_hash,
    uint32_t weight,
    const std::vector<double> &probability_x4,
    uint32_t seq_length,
    uint32_t seq_count);

// Предарссчет вероятностей для 4 букв: probability_x4
std::vector<double> calc_probability_x4(
    const std::vector<std::string> &sequences,
    bool search_complementary);

// Частоты для нуклеотидов (с учетом или без комплементарности)
std::vector<double> calc_frequencies_ratio(const std::vector<std::string> &sequences, bool complementary);

// Встречаемость нуклеотидов
std::vector<uint32_t> calc_frequencies(const std::vector<std::string> &sequences, bool complementary = false);

// найти такое p, после которого вероятность встретить по случайным причиныам меньше max_value (с помощью мат ожиданя)
double calc_prob_upper_bound_simple(double external_percentage, unsigned int positions);
std::vector<double> calc_probability_x4(const std::vector<double> &frequencies_ratio);