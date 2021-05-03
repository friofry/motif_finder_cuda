#include "gmps.h"

double gmps(int time, int seq_len, int n_str)
{
    return 2.562890625 * seq_len * n_str / double(time);
}


double gmps_external(double time, int number_of_hashes, unsigned int motifs_count)
{
    double ratio = double(motifs_count) * number_of_hashes / 1000000000.0;
    return ratio / double(time);
}

double gmps_external(double time, int seq_len, int seq_count, bool complementary, unsigned int motifs_count)
{
    unsigned int number_of_hashes = (seq_len - 7) * seq_count * (complementary ? 2 : 1);
    return gmps_external(time, number_of_hashes, motifs_count);
}