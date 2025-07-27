#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "baligner.hpp"


struct Anchor {
    uint query_start;
    uint ref_start;
};

std::vector<OpLen> merge_cigar_elements(const std::vector<OpLen>& elements) {
    if (elements.empty()) {
        return {};
    }
    std::vector<OpLen> merged_elements;
    merged_elements.push_back(elements[0]);
    for (size_t i = 1; i < elements.size(); ++i) {
        if (elements[i].op == merged_elements.back().op) {
            merged_elements.back().len += elements[i].len;
        } else {
            merged_elements.push_back(elements[i]);
        }
    }
    return merged_elements;
}

AlignmentResult piecewise_extension_alignment(
    const std::string& query,
    const std::string& reference,
    const std::vector<Anchor>& anchors,
    const int k,
    const int padding,
    const AlignmentScoring& scoring_params
) {
    AlignmentResult result;
    result.score = 0;
    std::vector<OpLen> temp_cigar_elements;

    const Anchor& first_anchor = anchors[0];
    if (first_anchor.query_start > 0 && first_anchor.ref_start > 0) {
        std::string query_part = query.substr(0, first_anchor.query_start);
        const size_t ref_start = std::max(0, static_cast<int>(first_anchor.ref_start) - (static_cast<int>(query_part.length()) + padding));
        std::string ref_part = reference.substr(ref_start, first_anchor.ref_start - ref_start);

        AlignmentResult pre_align = free_query_start_alignment(query_part, ref_part, scoring_params);
        result.score += pre_align.score;
        result.query_start = 0;
        result.ref_start = ref_start + pre_align.ref_start;
        temp_cigar_elements.insert(temp_cigar_elements.end(), pre_align.cigar.begin(), pre_align.cigar.end());
    } else {
        result.query_start = first_anchor.query_start;
        result.ref_start = first_anchor.ref_start;
    }

    result.score += k * scoring_params.match;
    temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)k});

    for (size_t i = 1; i < anchors.size(); ++i) {
        const Anchor& anchor = anchors[i];
        const Anchor& prev_anchor = anchors[i - 1];

        size_t curr_start_query = anchor.query_start;
        size_t curr_start_ref = anchor.ref_start;
        size_t prev_end_query = prev_anchor.query_start + k;
        size_t prev_end_ref = prev_anchor.ref_start + k;

        size_t ref_diff = curr_start_ref - prev_end_ref;
        size_t query_diff = curr_start_query - prev_end_query;

        if (ref_diff > 0 || query_diff > 0){
            std::string query_part = query.substr(prev_end_query, query_diff);
            std::string ref_part = reference.substr(prev_end_ref, ref_diff);

            AlignmentResult aligned = global_alignment(query_part, ref_part, scoring_params);
            result.score += aligned.score;
            temp_cigar_elements.insert(temp_cigar_elements.end(), aligned.cigar.begin(), aligned.cigar.end());

            result.score += k * scoring_params.match;
            temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)k});
        } else {
             if (ref_diff < query_diff) {
                const size_t inserted_part = ref_diff - query_diff;
                result.score += scoring_params.gap_open + (inserted_part - 1) * scoring_params.gap_extend;
                temp_cigar_elements.push_back({Operation::I, (uintptr_t)inserted_part});

                const size_t matching_part = k + ref_diff;
                result.score += matching_part * scoring_params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            } else if (ref_diff > query_diff) {
                const size_t deleted_part = query_diff - ref_diff;
                result.score += scoring_params.gap_open + (deleted_part - 1) * scoring_params.gap_extend;
                temp_cigar_elements.push_back({Operation::D, (uintptr_t)deleted_part});

                const size_t matching_part = k + query_diff;
                result.score += matching_part * scoring_params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            } else {
                const size_t matching_part = k + ref_diff;
                result.score += matching_part * scoring_params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            }
        }
    }

    const Anchor& last_anchor = anchors.back();
    const size_t last_anchor_end_query = last_anchor.query_start + k;
    const size_t last_anchor_end_ref = last_anchor.ref_start + k;
    if (last_anchor_end_query < query.length() && last_anchor_end_ref < reference.length()) {
        std::string query_part = query.substr(last_anchor_end_query);
        const size_t ref_part_end = std::min(reference.length(), last_anchor_end_ref + query_part.length() + padding);
        std::string ref_part = reference.substr(last_anchor_end_ref, ref_part_end - last_anchor_end_ref);

        AlignmentResult post_align = free_query_end_alignment(query_part, ref_part, scoring_params);
        result.score += post_align.score;
        result.query_end = last_anchor_end_query + post_align.query_end;
        result.ref_end = last_anchor_end_ref + post_align.ref_end;
        temp_cigar_elements.insert(temp_cigar_elements.end(), post_align.cigar.begin(), post_align.cigar.end());
    } else {
        result.query_end = last_anchor_end_query;
        result.ref_end = last_anchor_end_ref;
    }

    result.cigar = merge_cigar_elements(temp_cigar_elements);
    return result;
}

void visualize_alignment(const std::string& query, const std::string& reference, const AlignmentResult& result) {
    std::stringstream aligned_query_ss;
    std::stringstream aligned_ref_ss;
    std::stringstream match_string_ss;

    size_t ref_pos = 0;

    while (ref_pos < result.ref_start) {
        aligned_ref_ss << reference[ref_pos];
        aligned_query_ss << '.';
        match_string_ss << ' ';
        ref_pos++;
    }

    size_t query_pos = result.query_start;

    for (const auto& op_element : result.cigar) {
        uint count = op_element.len;
        Operation operation = op_element.op;

        switch (operation) {
            case Operation::M:
            case Operation::Eq:
                for (uint i = 0; i < count; ++i) {
                    aligned_ref_ss << reference[ref_pos++];
                    aligned_query_ss << query[query_pos++];
                    match_string_ss << '|';
                }
                break;
            case Operation::X:
                for (uint i = 0; i < count; ++i) {
                    aligned_ref_ss << reference[ref_pos++];
                    aligned_query_ss << query[query_pos++];
                    match_string_ss << 'X';
                }
                break;
            case Operation::I:
                for (uint i = 0; i < count; ++i) {
                    aligned_ref_ss << '-';
                    aligned_query_ss << query[query_pos++];
                    match_string_ss << ' ';
                }
                break;
            case Operation::D:
                for (uint i = 0; i < count; ++i) {
                    aligned_ref_ss << reference[ref_pos++];
                    aligned_query_ss << '-';
                    match_string_ss << ' ';
                }
                break;
            case Operation::Sentinel:
                break;
        }
    }

    while (ref_pos < reference.length()) {
        aligned_ref_ss << reference[ref_pos];
        aligned_query_ss << '.';
        match_string_ss << ' ';
        ref_pos++;
    }

    std::cout << "Ref    : " << aligned_ref_ss.str() << std::endl;
    std::cout << "         " << match_string_ss.str() << std::endl;
    std::cout << "Query  : " << aligned_query_ss.str() << std::endl;
    std::cout << std::endl;
}


int main() {
    const std::string reference = "AAAATCACTACCCGCATACGTTCCCC";
    const std::string query = "TCACTAACCGCTACGAT";
    const int k = 3;
    const int padding = 2;

    const std::vector<Anchor> anchors = {
        {2, 6},
        {7, 11},
        {12, 17}
    };

    AlignmentScoring default_scoring = {
        .match = 3,
        .mismatch = -1,
        .gap_open = -3,
        .gap_extend = -1
    };

    AlignmentResult result = piecewise_extension_alignment(query, reference, anchors, k, padding, default_scoring);

    std::cout << "Piecewise Extension Alignment Test\n" << std::endl;
    std::cout << "Query: " << query << std::endl;
    std::cout << "Ref:   " << reference << std::endl;
    std::cout << "Score: " << result.score << std::endl;
    std::cout << "Query Range: " << result.query_start << " - " << result.query_end << std::endl;
    std::cout << "Ref Range:   " << result.ref_start << " - " << result.ref_end << std::endl;
    std::cout << "CIGAR: " << result.to_cigar_string() << std::endl;
    std::cout << std::endl;

    visualize_alignment(query, reference, result);

    return 0;
}
