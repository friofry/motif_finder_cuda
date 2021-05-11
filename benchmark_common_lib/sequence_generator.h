#ifndef MOTIF_FINDER_SEQUENCE_GENERATOR_H
#define MOTIF_FINDER_SEQUENCE_GENERATOR_H

#include <string>
#include <vector>
#include <cstdint>

std::vector<std::string> generate_random_sequences(uint32_t  seq_length, uint32_t  seq_count);
std::vector<std::string> fetch_sequences_from_file(uint32_t  seq_length, uint32_t  seq_count, std::string seq_file);

#endif // MOTIF_FINDER_SEQUENCE_GENERATOR_H
