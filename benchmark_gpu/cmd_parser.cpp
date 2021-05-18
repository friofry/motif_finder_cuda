//
// Created by andrey on 03.03.2021.
//
#include "cmd_parser.h"
#include <config.h>
#include <cxxopts.hpp>
#include <gmps.h>
#include <iostream>
#include <json.hpp>

using namespace std;
using json = nlohmann::json;

template<typename T>
T parse_param(const cxxopts::ParseResult &result, const char *name, const T &default_val)
{
    if (result.count(name) == 0) {
        return default_val;
    }
    return result[name].as<T>();
}

template<typename T>
void print_vec(const string &name, const vector<T> &v)
{
    cout << name << " ";
    for (const auto x: v) {
        cout << v << " ";
    }
    cout << endl;
}

vector<int> parse_int_vec(const cxxopts::ParseResult &result, const char *name, const vector<int> &default_val)
{
    auto res = parse_param<vector<int>>(result, name, default_val);
    //print_vec<vector<int>>(name, res);
    return res;
}
vector<string> parse_string_vec(const cxxopts::ParseResult &result, const char *name, const vector<string> &default_val)
{
    auto res = parse_param<vector<string>>(result, name, default_val);
    //print_vec<vector<string>>(name, res);
    return res;
}
string parse_string_param(const cxxopts::ParseResult &result, const char *name, const string &default_val)
{
    return parse_param<string>(result, name, default_val);
}
TestParams parse_gpu_input_args(int argc, char **argv)
{
    cxxopts::Options options("Test Benchmark", "Test gpu external and internal algorithms performance");
    options.add_options()
            ("l,lengths", "Sequence lengths", cxxopts::value<vector<int>>())
            ("c,counts", "Sequence counts", cxxopts::value<vector<int>>())
            ("g,gpu_counts", "GPU counts", cxxopts::value<vector<int>>())
            ("o,complementary", "Complementary values", cxxopts::value<vector<int>>())
            ("a,algorithms", "internal|external", cxxopts::value<vector<string>>())
            ("r,external_min_presence", "Maximum motif presence by chance (%)", cxxopts::value<vector<int>>())
            ("m,external_stat_model", "bernulli|ideal|markov1|markov2|markov3", cxxopts::value<vector<string>>())
            ("f,sequencesfile", "File with sequences", cxxopts::value<string>())
            ("h,threads_per_block", "Threads per block", cxxopts::value<vector<int>>())
            ("u,chunk_sizes", "Motif chunk sizes", cxxopts::value<vector<int>>())
            ("i,unified_memory", "Use unified memory API", cxxopts::value<vector<int>>())
            ;

    auto r = options.parse(argc, argv);

    TestParams res;
    res.seq_counts = parse_int_vec(r, "counts", {10});
    res.seq_lengths = parse_int_vec(r, "lengths", {10});
    res.gpu_counts = parse_int_vec(r, "gpu_counts", {1});
    res.complementary = parse_int_vec(r, "complementary", {0});
    res.algorithms = parse_string_vec(r, "algorithms", {"internal"});
    res.sequences_file = parse_string_param(r, "sequencesfile", "");
    res.threads_per_block = parse_int_vec(r, "threads_per_block", {THREADS_PER_BLOCK});
    res.chunk_sizes = parse_int_vec(r, "chunk_sizes", {MOT_PER_CHUNK});
    res.external_min_presence = parse_int_vec(r, "external_min_presence", {100});
    res.external_stat_model = parse_string_vec(r, "external_stat_model", {"ideal"});
    res.unified_memory = parse_int_vec(r, "unified_memory", {1});

    return res;
}

std::string gpu_result_json(const GpuRunParams &params, float result, unsigned int motifs_count)
{
    json o;
    o["count"] = params.count;
    o["length"] = params.length;
    o["gpus"] = params.gpus;
    o["complementary"] = params.complementary;
    o["algorithm"] = params.algorithm;
    o["threads_per_block"] = params.threads_per_block;
    o["motif_chunk_size"] = params.motif_chunk_size;
    o["external_min_presence"] = params.external_min_presence;
    o["external_stat_model"] = params.external_stat_model;
    o["unified_memory"] = params.unified_memory;
    o["result"] = result;
    o["motifs_count"] = motifs_count;
    o["gmps"] = gmps_external(result, params.length, params.count, params.complementary, motifs_count);
    o["gmps_total"] = gmps_external(result, params.length, params.count, params.complementary, TOTAL_MOT);
    auto res = o.dump();
    //printf("%s\n",res.c_str());
    return res;
}

string combine_results_json(const std::vector<std::string> &results)
{
    json ret;
    for (const auto &r: results) {
        auto j = json::parse(r);
        ret.push_back(j);
    }
    return ret.dump();
}
