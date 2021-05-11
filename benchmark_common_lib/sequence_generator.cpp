#include "sequence_generator.h"

#include <letter_conversions.h>
#include <fst_reader.h>

using namespace std;

std::vector<std::string> generate_random_sequences(uint32_t seq_length, uint32_t seq_count)
{
    vector<string> sequences;

    for (uint32_t  i = 0; i < seq_count; i++) {
        string sequence;
        for (uint32_t  j = 0; j < seq_length; j++) {
            sequence += to_symbol(rand() % 4);
        }
        sequences.push_back(sequence);
    }
    return sequences;
}

std::vector<std::string> fetch_sequences_from_file(uint32_t seq_length, uint32_t seq_count, std::string seq_file)
{
    if (seq_file.empty()) {
        return generate_random_sequences(seq_length, seq_count);
    }

    vector<string> file_sequences;

    try {
        file_sequences = read_fasta(seq_file.c_str());
    } catch (const invalid_argument &ex) {
    }

    vector<string> sequences;

    if (file_sequences.size() < seq_count) {
        return generate_random_sequences(seq_length, seq_count);
    }

    for (uint32_t  i = 0; i < seq_count; i++) {
        string sequence;
        const auto &file_sequence = file_sequences[i];
        if (file_sequence.length() < seq_length) {
            return generate_random_sequences(seq_length, seq_count);
        }
        for (uint32_t  j = 0; j < seq_length; j++) {
            sequence += file_sequence[j];
        }
        sequences.push_back(sequence);
    }
    return sequences;
}

