#pragma once

double gmps(int time, int seq_len, int n_str);
double gmps_external(double time, int number_of_hashes, unsigned int motifs_count);
double gmps_external(double time, int seq_len, int seq_count, bool complementary, unsigned int motifs_count);
