#include <vector>
#include <cstdint>

#include <InternalCpuAlgorithm.h>
#include <ExternalCpuAlgorithm.h>

#include <fst_reader.h>
#include <sequences_to_hashes.h>

using namespace std;

void run()
{
    const auto sequences = read_fasta("test.fst");
    SequenceHashes hashesInfo = sequencesToHashes(sequences, false);
    std::vector<uint32_t> hashes { 1, 2, 3};
    std::vector<uint16_t> res;
    ExternalCpuAlgorithm external(hashes, hashesInfo);
    external.getAllOccurrences(res);
//    InternalCpuAlgorithm internal(hashesInfo);
//    internal.getAllOccurrences(res);
}

int main(int argc, char **argv)
{
    argc = 0;
    argv = 0;
    run();
    return 0;
}
