#include "fst_reader.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <regex>

using namespace std;

namespace {

bool match_fasta_sequence(const string &s)
{
    regex rgx("[atgcATGC]+");
    return regex_match(s, rgx);
}

bool match_fasta_header(const string &s)
{
    regex rgx("\\s*>.*");
    return regex_match(s, rgx);
}

std::string trim(const std::string &str)
{
    const char *ws = " \t\n\r\f\v";
    string result = str;
    result.erase(0, str.find_first_not_of(ws));
    result.erase(result.find_last_not_of(ws)+1);
    return result;
}

string find_sequence_id(const string &header)
{
    auto pos = header.find(">");
    if (pos == std::string::npos) {
        return "";
    }
    return trim(header.substr(pos+1));
}

bool match_whitespace(const string &s)
{
    regex rgx("\\s*");
    return regex_match(s, rgx);
}

} // namespace

vector<string> read_fasta(const char *filename)
{
    vector<string> result;
    if (strlen(filename) == 0) {
        return result;
    }

    ifstream input(filename);

    if (input.eof() || input.bad() || !input.good()) {
        throw invalid_argument(std::string("Cannot open input file: ") + filename);
        return result;
    }

    string line;
    string id;
    string current_sequence;
    int line_number = 0;

    while (getline(input, line)) {
        line_number++;
        line = trim(line);

        if (line.size() == 0 || match_whitespace(line)) {
            continue;
        } else if (match_fasta_header(line)) {
            if (!current_sequence.empty()) {
                result.push_back(current_sequence);
                current_sequence.clear();
            }
            id = find_sequence_id(line);
        } else {
            if (!match_fasta_sequence(line)) {
                throw invalid_argument("Invalid character in line:" + to_string(line_number) + ":\n" + line);
            };
            current_sequence += trim(line);
        }
    }

    if (!current_sequence.empty()) {
        result.push_back(current_sequence);
    }
    return result;
}
