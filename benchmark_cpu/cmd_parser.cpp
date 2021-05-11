#include "cmd_parser.h"
#include <config.h>
#include <cxxopts.hpp>
#include <gmps.h>
#include <iostream>
#include <json.hpp>

#include <parse_helper.h>

using namespace std;
using json = nlohmann::json;

CpuTestParams parse_input_arguments(int argc, char **argv)
{
    cxxopts::Options options("Test Benchmark", "Test gpu external and internal algorithms performance");
    options.add_options()
            ("l,lengths", "Sequence lengths", cxxopts::value<vector<int>>())
            ("c,counts", "Sequence counts", cxxopts::value<vector<int>>())
            ("j,cpu_threads", "CPU threads", cxxopts::value<vector<int>>())
            ("o,complementary", "Complementary values", cxxopts::value<vector<int>>())
            ("a,algorithms", "internal|external", cxxopts::value<vector<string>>())
            ("r,external_min_presence", "Maximum motif presence by chance (%)", cxxopts::value<vector<int>>())
            ("m,external_stat_model", "bernulli|ideal|markov1|markov2|markov3", cxxopts::value<vector<string>>())
            ("f,sequencesfile", "File with sequences", cxxopts::value<string>())
            ;

    auto r = options.parse(argc, argv);

    CpuTestParams res;
    res.seq_counts = parse_int_vec(r, "counts", {10});
    res.seq_lengths = parse_int_vec(r, "lengths", {10});
    res.cpu_threads = parse_int_vec(r, "cpu_threads", {1});
    res.complementary = parse_int_vec(r, "complementary", {0});
    res.algorithms = parse_string_vec(r, "algorithms", {"internal"});
    res.sequences_file = parse_string_param(r, "sequencesfile", "");
    res.external_min_presence = parse_int_vec(r, "external_min_presence", {100});
    res.external_stat_model = parse_string_vec(r, "external_stat_model", {"ideal"});

    return res;
}

std::string make_result_json(const CpuRunParams &params, float result, unsigned int motifs_count)
{
    json o;
    o["count"] = params.count;
    o["length"] = params.length;
    o["cpu_threads"] = params.cpu_threads;
    o["complementary"] = params.complementary;
    o["algorithm"] = params.algorithm;
    o["external_min_presence"] = params.external_min_presence;
    o["external_stat_model"] = params.external_stat_model;
    o["result"] = result;
    o["motifs_count"] = motifs_count;
    o["gmps"] = gmps_external(result, params.length, params.count, params.complementary, motifs_count);
    o["gmps_total"] = gmps_external(result, params.length, params.count, params.complementary, TOTAL_MOT);
    auto res = o.dump();
    //printf("%s\n",res.c_str());
    return res;
}

