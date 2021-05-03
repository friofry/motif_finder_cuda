#include <sstream>
#include <vector>
#include <cstdint>

#include <motif_finder_cpu_internal.h>
#include <fst_reader.h>

using namespace std;

void run()
{
    CpuInternalParams params;
    params.sequences = read_fasta("test.fst");
    params.search_complementary = false;

    std::vector<uint16_t> res;
    find_motifs_internal_cpu(res, params);
}

int main(int argc, char **argv)
{
    argc = 0;
    argv = 0;
    run();
    return 0;
}
