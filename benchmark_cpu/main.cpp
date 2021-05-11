#include <vector>
#include <cstdint>

#include <internal_cpu_algorithm.h>
#include <external_cpu_algorithm.h>

#include <fst_reader.h>
#include <sequences_to_hashes.h>

using namespace std;

void run()
{
    const auto sequences = read_fasta("test.fst");
    SequenceHashes sequence_hashes = sequences_to_hashes(sequences, false);
    std::vector<uint16_t> occurrences;
    internal_cpu_algorithm(sequence_hashes, occurrences);
}

int main(int argc, char **argv)
{
    argc = 0;
    argv = 0;
    run();
    return 0;
}
