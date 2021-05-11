#ifndef MOTIF_FINDER_CMD_PARSER_H
#define MOTIF_FINDER_CMD_PARSER_H
#include <vector>
#include <string>

struct TestParams {
    std::vector<int> seq_counts;
    std::vector<int> seq_lengths;
    std::vector<int> gpu_counts;
    std::vector<int> complementary;
    std::vector<std::string> algorithms;
    std::string sequences_file;
    std::vector<int> threads_per_block;
    std::vector<int> chunk_sizes;
    std::vector<std::string> external_stat_model;
    std::vector<int> external_min_presence;
    std::vector<int> unified_memory;
};

struct RunParams {
    int count {0};
    int length{0};
    int gpus{0};
    bool complementary{false};
    std::string algorithm;
    int threads_per_block{0};
    int motif_chunk_size{0};
    std::string external_stat_model;
    int external_min_presence{0};
    bool unified_memory{true};
};

TestParams parse_gpu_input_args(int argc, char **argv);
std::string gpu_result_json(const RunParams &params, float result, unsigned int motifs_count);
std::string combine_results_json(const std::vector<std::string> &results);

#endif //MOTIF_FINDER_CMD_PARSER_H
