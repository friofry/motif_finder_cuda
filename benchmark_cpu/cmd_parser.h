#ifndef CPU_MOTIF_FINDER_CMD_PARSER_H
#define CPU_MOTIF_FINDER_CMD_PARSER_H
#include <vector>
#include <string>

struct CpuTestParams {
    std::vector<int> seq_counts;
    std::vector<int> seq_lengths;
    std::vector<int> cpu_threads;
    std::vector<int> complementary;
    std::vector<std::string> algorithms;
    std::string sequences_file;
    std::vector<std::string> external_stat_model;
    std::vector<int> external_min_presence;
};

struct CpuRunParams {
    int count {0};
    int length {0};
    int cpu_threads {1};
    bool complementary{false};
    std::string algorithm;
    std::string external_stat_model;
    int external_min_presence {0};
};

CpuTestParams parse_input_arguments(int argc, char **argv);
std::string make_result_json(const CpuRunParams &params, float result, unsigned int motifs_count);

#endif //CPU_MOTIF_FINDER_CMD_PARSER_H
