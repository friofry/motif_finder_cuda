#pragma once

#include <config.h>
#include <cstdint>
#include <stat_model.h>
#include <string>
#include <vector>

// Получить все индексы мотивов, для которых вероятность встретить по случайным причинам меньше заданной
// С учетом частот во входных последовательностях
void generate_external_indexes_cpu(double external_percentage,
                                   const std::vector<std::string> &sequences,
                                   bool complementary,
                                   std::vector<uint32_t> &result);

// Получить все индексы мотивов, для которых вероятность встретить по случайным причинам меньше заданной
void generate_external_indexes_cpu(double max_motif_prob_simple,
                                   const std::vector<double> &probability_x4,
                                   std::vector<uint32_t> &result);

// Сгенерировать список мотивов с помощью заданной StatModel
void generate_external_hashes_cpu(double max_motif_prob_by_chance,
                                  const StatModel &stat_model,
                                  bool complementary,
                                  std::vector<uint32_t> &result_hashes);

// Отфильтровать список мотивов с помощью контрастной выборки
void filter_hashes_by_contrast(double max_score,
                               const StatModel &stat_model,
                               const std::vector<uint16_t> &weights,
                               std::vector<uint32_t> &motifs_hashes);

void generate_all_hashes(std::vector<uint32_t> &motif_hashes, uint32_t count = TOTAL_MOT);