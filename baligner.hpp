#ifndef BLOCK_ALIGNER_WRAPPER_H
#define BLOCK_ALIGNER_WRAPPER_H
#include <string>
#include <vector>
#include <cstddef>
#include "block_aligner.h"

struct AlignmentScoring {
    int8_t match;
    int8_t mismatch;
    int8_t gap_open;
    int8_t gap_extend;
};

struct AlignmentResult {
    int score;
    size_t query_start;
    size_t query_end;
    size_t ref_start;
    size_t ref_end;
    std::vector<OpLen> cigar;
    std::string to_cigar_string() const;
};

AlignmentResult global_alignment(const std::string& query, const std::string& ref, const AlignmentScoring& scoring_params);
AlignmentResult free_query_end_alignment(const std::string& query, const std::string& ref, const AlignmentScoring& scoring_params);
AlignmentResult free_query_start_alignment(const std::string& query, const std::string& ref, const AlignmentScoring& scoring_params);


#endif
