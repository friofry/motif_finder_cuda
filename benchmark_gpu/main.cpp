#include <sstream>
#include <vector>
#include <cstdint>

// cpu_lib
#include <external_gpu_algorithm.h>
#include <internal_gpu_algorithm.h>
#include <gpu_cuda_params.h>

// common
#include <fst_reader.h>
#include <gmps.h>
#include <timer.h>

// benchmark_common
#include <benchmark_common.h>
#include <parse_helper.h>
#include <sequence_generator.h>

// conversion
#include <sequences_to_hashes.h>

// stat
#include <external_generator_cpu.h>
#include <stat_model.h>

#include "cmd_parser.h"

using namespace std;

std::string run(const vector<uint32_t> &motif_hashes, const SequenceHashes &sequence_hashes, const GpuRunParams &params)
{
    Timer t;
    t.silence();

    GpuCudaParams cuda_params;
    cuda_params.gpu_count = params.gpus;
    cuda_params.unified_memory = params.unified_memory;
    cuda_params.threads_per_block = params.threads_per_block;
    cuda_params.motif_range_size = params.motif_chunk_size;

    if (params.algorithm == "internal") {
        std::vector<uint16_t> occurrences;
        internal_gpu_algorithm(sequence_hashes, occurrences, cuda_params);
    }

    if (params.algorithm == "external") {
        std::vector<uint16_t> occurrences;
        external_gpu_algorithm(motif_hashes,sequence_hashes, occurrences, cuda_params);
    }

    float result = t.elapsed();
    std::stringstream ss;
    ss << params.algorithm.c_str() << "\t"
       << params.complementary << "\t"
       << params.gpus << "\t"
       << params.count << "\t"
       << params.length << "\t"
       << params.threads_per_block << "\t"
       << params.motif_chunk_size << "\t"
       << params.external_stat_model << "\t"
       << params.external_min_presence << "\t"
       << result;
    auto res = ss.str();

    auto motif_count = motif_hashes.size();
    if (params.algorithm == "internal") {
        motif_count = TOTAL_MOT;
    }

    FILE * f = fopen ("_results.txt","a+");
    fprintf(f,"%s\n", res.c_str());
    fclose(f);

    double gmps = gmps_external(result, params.length, params.count, params.complementary, motif_count);
    printf("%s %f\n", res.c_str(), gmps);

    return gpu_result_json(params, result, motif_count);
}

int main(int argc, char **argv)
{
    auto params = parse_gpu_input_args(argc, argv);

    vector<string> whole_res;
    bool only_internal = params.algorithms.size() == 1 && params.algorithms[0] == "internal";

    for (int length : params.seq_lengths) {
        for (int count : params.seq_counts) {
            const auto &sequences = fetch_sequences_from_file(length, count, params.sequences_file);
            stringstream ss;
            ss << length << "_" << count << ".fst";
            save_to_file(sequences, ss.str().c_str());

            for (int complementary : params.complementary) {
                const auto &sequence_hashes = sequences_to_hashes(sequences, complementary);
                for (int external_presence : params.external_min_presence) {
                    for (auto external_model : params.external_stat_model) {
                        vector<uint32_t> motif_hashes;
                        StatModelPtr stat_model = create_stat_model(external_model, sequences, complementary);

                        // Если не получилось считать хэши, рассчитать заново
                        if (!only_internal) {
                            generate_external_hashes_cpu(
                                external_presence * 0.01, *stat_model, complementary, motif_hashes);
                        }

                        for (auto algorithm : params.algorithms) {
                            for (int gpu_count : params.gpu_counts) {
                                for (int threads_per_block : params.threads_per_block) {
                                    for (int chunk_size : params.chunk_sizes) {
                                        for (bool unified_mem : params.unified_memory) {
                                            GpuRunParams cur_params;
                                            cur_params.count = count;
                                            cur_params.length = length;
                                            cur_params.gpus = gpu_count;
                                            cur_params.complementary = complementary;
                                            cur_params.algorithm = algorithm;
                                            cur_params.threads_per_block = threads_per_block;
                                            cur_params.motif_chunk_size = chunk_size;
                                            cur_params.external_min_presence = external_presence;
                                            cur_params.external_stat_model = external_model;
                                            cur_params.unified_memory = unified_mem;

                                            auto res = run(motif_hashes, sequence_hashes, cur_params);
                                            whole_res.push_back(res);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    auto result_json = combine_results_json(whole_res);
    FILE *f = fopen("results.json", "w+");
    fprintf(f, "%s\n", result_json.c_str());
    fclose(f);

    return 0;
}