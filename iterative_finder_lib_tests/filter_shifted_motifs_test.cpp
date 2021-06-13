//
// Created by andrey on 13.06.2021.
//
#include <gtest/gtest.h>
#include <hash_conversions.h>
#include "filter_shifted_motifs.h"

bool matches_to_shift(const char *motif, const char *query, bool complementary) {
    auto motif_hash = string_to_hash(motif);
    auto query_hash = string_to_hash(query);
    return hash_matches(query_hash, {motif_hash}, complementary);
}

TEST(FilterShiftedTests, FiltersRightShift) {
    EXPECT_TRUE(matches_to_shift("AATGGGGG", "NAATGGGG", false));
}

TEST(FilterShiftedTests, NotFiltersRightShift) {
    EXPECT_FALSE(matches_to_shift("AATGGGGG", "WAATGGGG", false));
}

TEST(FilterShiftedTests, FilterLeftShift) {
    EXPECT_TRUE(matches_to_shift("AATGGGGG", "ATGGGGGN", false));
}

TEST(FilterShiftedTests, NotFilterLeftShift) {
    EXPECT_FALSE(matches_to_shift("AATGGGGG", "ATGGGGGW", false));
}

TEST(FilterShiftedTests, FiltersRightShiftComplementary) {
    EXPECT_TRUE(matches_to_shift("CCCCCATT", "NAATGGGG", true));
}

TEST(FilterShiftedTests, NotFiltersRightShiftComplementary) {
    EXPECT_FALSE(matches_to_shift("CCCCCATT", "WAATGGGG", true));
}

TEST(FilterShiftedTests, FilterLeftShiftComplementary) {
    EXPECT_TRUE(matches_to_shift("CCCCCATT", "ATGGGGGN", true));
}

TEST(FilterShiftedTests, NotFilterLeftShiftComplementary) {
    EXPECT_FALSE(matches_to_shift("CCCCCATT", "ATGGGGGB", true));
}

TEST(FilterShiftedTests, NotFiltersEqual) {
    EXPECT_FALSE(matches_to_shift("CCCCCATT", "CCCCCATT", false));
    EXPECT_FALSE(matches_to_shift("CCCCCATT", "CCCCCATT", true));
}