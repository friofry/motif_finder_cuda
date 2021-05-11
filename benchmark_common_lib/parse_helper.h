
#ifndef MOTIF_FINDER_PARSE_HELPER_H
#define MOTIF_FINDER_PARSE_HELPER_H

#include <cxxopts.hpp>
#include <json.hpp>
#include <vector>
#include <string>

std::vector<int> parse_int_vec(const cxxopts::ParseResult &result, const char *name, const std::vector<int> &default_val);
std::vector<std::string> parse_string_vec(const cxxopts::ParseResult &result, const char *name, const std::vector<std::string> &default_val);
std::string parse_string_param(const cxxopts::ParseResult &result, const char *name, const std::string &default_val);

std::string combine_results_json(const std::vector<std::string> &results);

#endif //MOTIF_FINDER_PARSE_HELPER_H
