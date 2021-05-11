#include <cxxopts.hpp>
#include <json.hpp>
#include <iostream>

using namespace std;
using json = nlohmann::json;

namespace {
    template<typename T>
    T parse_param(const cxxopts::ParseResult &result, const char *name, const T &default_val) {
        if (result.count(name) == 0) {
            return default_val;
        }
        return result[name].as<T>();
    }

    template<typename T>
    void print_vec(const string &name, const vector<T> &v) {
        cout << name << " ";
        for (const auto x: v) {
            cout << v << " ";
        }
        cout << endl;
    }
} // namespace

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

string combine_results_json(const std::vector<std::string> &results)
{
    json ret;
    for (const auto &r: results) {
        auto j = json::parse(r);
        ret.push_back(j);
    }
    return ret.dump();
}
