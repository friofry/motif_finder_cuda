#include "markov_stat_model.h"

#include <unordered_map>
#include <functional>

#include <config.h>
#include <letter_conversions.h>
#include <fst_reader.h>
#include <timer.h>

#include "probability.h"

using namespace std;

namespace {

const double _default_prob = 0.000000001;

uint32_t number_of_hashes(int kmer_length)
{
    return 1 << (4 * kmer_length);
}

const vector<vector<int>> &iupac_codes_by_nucl()
{
    static vector<vector<int>> result;

    if (result.empty()) {
        result.resize(4);
        for (int i = 0; i < 4; i++) {
            uint32_t hash = num_to_hash(i);
            vector<int> iupac_codes;
            for (int j = 1; j < 16; j++) {
                if (j & hash) {
                    iupac_codes.push_back(j);
                }
            }
            result[i] = iupac_codes;
        }
    }

    return result;
}

bool has_key(const unordered_map<uint32_t, int> &counters, uint32_t key)
{
    return counters.find(key) != counters.end();
}

void increment(unordered_map<uint32_t, int> &counters, uint32_t key)
{
    if (!has_key(counters, key)) {
        counters[key] = 0;
    }
    counters[key]++;
}

uint32_t make_hash(const uint8_t *start, int count)
{
    uint32_t hash = 0;
    const uint8_t *cur = start;
    for (int i = 0; i < count; i++) {
        hash <<= 4;
        hash |= num_to_hash(*cur);
        cur++;
    }
    return hash;
}

void increment_nucl_hashes(const SequenceNums &sequence_nums, unordered_map<uint32_t, int> &hash_counters, uint32_t kmer_length)
{
    int offset = 0;
    for (uint32_t i = 0; i < sequence_nums.count; i++) {
        int length = sequence_nums.lengths[i];
        for (uint32_t j = 0; j < length - kmer_length + 1; j++) {
            const uint8_t *start = &(sequence_nums.nums[offset + j]);
            uint32_t key_hash = make_hash(start, kmer_length);
            increment(hash_counters, key_hash);
        }
        offset += length;
    }
}

void iterate_slots(int positions, int max_in_slot, function<bool(const vector<int> &)> action)
{
    vector<int> slots(positions, 0);
    int index = 0;
    while (true) {
        bool ok = action(slots);
        if (!ok) {
            return;
        }
        slots[0]++;
        while (slots[index] == max_in_slot) {
            if (index == positions - 1) {
                return;
            }
            slots[index++] = 0;
            slots[index]++;
        }
        index = 0;
    }
}

uint32_t slots_to_hash(const vector<int> &slots, uint32_t nucl_hash)
{
    const auto &codes_by_nucl = iupac_codes_by_nucl();
    uint32_t cur_hash = 0;

    for (uint32_t i = 0; i < slots.size(); i++) {
        auto cur_num = nucl_hash_to_num(nucl_hash & 0xF);

        auto cur_iupac_hash = codes_by_nucl[cur_num][slots[i]];
        cur_hash |= cur_iupac_hash << (4 * i);

        nucl_hash >>= 4;
    }
    return cur_hash;
}

vector<uint32_t> calc_iupac_hashes_counts(const unordered_map<uint32_t, int> &nucl_hash_counters, int kmer_length)
{
    uint32_t hashes_count = number_of_hashes(kmer_length);
    vector<uint32_t> iupac_hashes_counter(hashes_count, 0);
    int iupacs_per_nucl = iupac_codes_by_nucl()[0].size();

    for (const auto& kv: nucl_hash_counters) {
        auto &nucl_hash = kv.first;
        auto &count = kv.second;

        auto update_all_iupacs = [&iupac_hashes_counter, count, nucl_hash](const vector<int> &slots) -> bool
        {
            uint32_t cur_hash = slots_to_hash(slots, nucl_hash);
            if (iupac_hashes_counter.size() <= cur_hash) {
                printf("error hash in calc_iupac_hashes_counts. size %lu, hash %d\n",  iupac_hashes_counter.size(), cur_hash);
                return false;
            }
            iupac_hashes_counter[cur_hash] += count;
            return true;
        };

        // Перебор для каждой позиции всех вариантов iupac кодов.
        iterate_slots(kmer_length, iupacs_per_nucl, update_all_iupacs);
    }
    return iupac_hashes_counter;
}


int kmers_in_sequences_by_length(const SequenceNums &sequence_nums, int kmer_length)
{
    int result = 0;
    for (uint32_t i = 0; i < sequence_nums.count; i++) {
        uint32_t length = sequence_nums.lengths[i];
        result += length - kmer_length + 1;
    }
    return result;
}

uint32_t extract_kmer_hash(uint32_t hash, int pos, int length)
{
    // Сделать маску
    uint32_t mask = 0xFFFFFFFF;
    int shift = (MOTIV_LEN - length) * 4;
    mask >>= shift;

    // Выделить хэш и сдвинуть в начало
    int shift_pos = shift - pos * 4;
    mask <<= shift_pos;
    return (mask & hash) >> shift_pos;
}

} // namespace


MarkovStatModel::MarkovStatModel(const std::vector<std::string> &sequences, bool complementary, int level)
    : StatModel(sequences, complementary)
    , _level(level)
    , _kmer_length(_level + 1)
{
    Timer t("MarkovStatModel");
    t.silence();
    if (level < 1 || level > 4) {
        throw invalid_argument("Supported levels for MarkovStatModel [1-4]");
    }

    // === посчитать вероятности
    uint32_t hashes_count = number_of_hashes(_kmer_length);
    _probabilities.resize(hashes_count, _default_prob);

    calc_probabilities(_kmer_length);
    calc_probabilities(_kmer_length-1);

    // === предарссчитанные данные
    _symbols_in_key1 = (MOTIV_LEN - _kmer_length) / 2 + 1;
    _symbols_in_key2 = MOTIV_LEN - _kmer_length - _symbols_in_key1 + 1;
    _key_size1 = _symbols_in_key1 + _kmer_length - 1;
    _key_size2 = _symbols_in_key2 + _kmer_length - 1;

    precalc_probabilities(_symbols_in_key1, _precalc_prefix_probabilities, true);
    precalc_probabilities(_symbols_in_key2, _precalc_suffix_probabilities, false);
}

uint32_t slots_to_hash(const vector<int> &slots)
{
    uint32_t hash = 0;

    for (uint32_t i = 0; i < slots.size(); i++) {
        hash |= (1 + slots[i]) << (i * 4);
    }
    return hash;
}

void MarkovStatModel::precalc_probabilities(int symbols_in_key, std::vector<double> &result_probabilites, bool prefix)
{
    int key_size = symbols_in_key + _kmer_length - 1;

    uint32_t keys_count = number_of_hashes(key_size);
    if (result_probabilites.size() < keys_count) {
        result_probabilites.resize(keys_count, 0.0);
    }

    auto precalc_key = [this, symbols_in_key, key_size, prefix, &result_probabilites](const vector<int> &slots) -> bool
    {
        uint32_t key_hash = slots_to_hash(slots);

        double probability = 1.0;
        int offset = MOTIV_LEN - key_size;

        if (prefix) {
            uint32_t cur_hash = extract_kmer_hash(key_hash, offset, _kmer_length);
            probability *= _probabilities[cur_hash];
        }

        int start = prefix ? 1 : 0;

        for (int i = start; i < symbols_in_key; i++) {
            uint32_t cur_hash = extract_kmer_hash(key_hash, offset + i, _kmer_length);
            probability *= _probabilities[cur_hash];

            uint32_t cur_hash_1 = extract_kmer_hash(key_hash, offset + i, _kmer_length-1);
            probability /= _probabilities[cur_hash_1];
        }

        result_probabilites[key_hash] = probability;

        return true;
    };

    // Перебор для каждой позиции всех вариантов iupac кодов.
    iterate_slots(key_size, ALPH_SIZE, precalc_key);
}

void MarkovStatModel::calc_probabilities(uint32_t kmer_length)
{
    // Встречаемость нуклеотидов k-меров из входных последовательностях
    unordered_map<uint32_t, int> nucl_hash_counters;
    increment_nucl_hashes(_sequence_nums, nucl_hash_counters, kmer_length);

    // Встречаемость iupac k-меров через встречаемость k-меров нуклеотидных
    vector<uint32_t> iupac_hashes_counters = calc_iupac_hashes_counts(nucl_hash_counters, kmer_length);

    // Вероятность iupac k-меров
    double total_kmers_count = static_cast<double>(kmers_in_sequences_by_length(_sequence_nums, kmer_length));

    for (uint32_t hash = 0; hash < iupac_hashes_counters.size(); hash++) {
        if (iupac_hashes_counters[hash] > 0) {
            if (hash >= _probabilities.size()) {
                printf("error hash in calc_probabilities\n");
            }
            _probabilities[hash] = iupac_hashes_counters[hash] / total_kmers_count;
        }
    }
}

double MarkovStatModel::motif_probability(uint32_t hash) const
{
    double probability = 1.0;

    for (int i = 0; i < MOTIV_LEN - _kmer_length + 1; i++) {
        uint32_t cur_hash = extract_kmer_hash(hash, i, _kmer_length);
        probability *= _probabilities[cur_hash];
    }

    double divisor = 1.0;
    for (int i = 1; i < MOTIV_LEN - _kmer_length + 1; i++) {
        uint32_t cur_hash_1 = extract_kmer_hash(hash, i, _kmer_length - 1);
        divisor *= _probabilities[cur_hash_1];
    }
    probability /= divisor;

    return probability;
}

double MarkovStatModel::motif_probability_x4(uint32_t hash) const
{
    double probability = 1.0;

    uint32_t key1 = extract_kmer_hash(hash, 0, _key_size1);
    uint32_t key2 = extract_kmer_hash(hash, _symbols_in_key1, _key_size2);

    probability *= _precalc_prefix_probabilities[key1];
    probability *= _precalc_suffix_probabilities[key2];
    return probability;
}