#include <algorithm>
#include <functional>
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

        if (pre_align.score == 0) {
            result.query_start = first_anchor.query_start;
            result.ref_start = first_anchor.ref_start;
        } else {
            result.score += pre_align.score;
            result.query_start = pre_align.query_start;
            result.ref_start = ref_start + pre_align.ref_start;
            temp_cigar_elements.insert(temp_cigar_elements.end(), pre_align.cigar.begin(), pre_align.cigar.end());
        }
    } else {
        result.query_start = first_anchor.query_start;
        result.ref_start = first_anchor.ref_start;
    }

    result.score += k * scoring_params.match;
    temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)k});

    for (size_t i = 1; i < anchors.size(); ++i) {
        const Anchor& anchor = anchors[i];
        const Anchor& prev_anchor = anchors[i - 1];

        int curr_start_query = anchor.query_start;
        int curr_start_ref = anchor.ref_start;
        int prev_end_query = prev_anchor.query_start + k;
        int prev_end_ref = prev_anchor.ref_start + k;

        int ref_diff = curr_start_ref - prev_end_ref;
        int query_diff = curr_start_query - prev_end_query;

        if (ref_diff > 0 && query_diff > 0){
            std::string query_part = query.substr(prev_end_query, query_diff);
            std::string ref_part = reference.substr(prev_end_ref, ref_diff);

            AlignmentResult aligned = global_alignment(query_part, ref_part, scoring_params);
            result.score += aligned.score;
            temp_cigar_elements.insert(temp_cigar_elements.end(), aligned.cigar.begin(), aligned.cigar.end());

            result.score += k * scoring_params.match;
            temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)k});
        } else {
             if (ref_diff < query_diff) {
                const size_t inserted_part = -ref_diff + query_diff;
                result.score += scoring_params.gap_open + (inserted_part - 1) * scoring_params.gap_extend;
                temp_cigar_elements.push_back({Operation::I, (uintptr_t)inserted_part});

                const size_t matching_part = k + ref_diff;
                result.score += matching_part * scoring_params.match;
                temp_cigar_elements.push_back({Operation::Eq, (uintptr_t)matching_part});
            } else if (ref_diff > query_diff) {
                const size_t deleted_part = -query_diff + ref_diff;
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

        if (post_align.score == 0) {
            result.query_end = last_anchor_end_query;
            result.ref_end = last_anchor_end_ref;
        } else {
            result.score += post_align.score;
            result.query_end = last_anchor_end_query + post_align.query_end;
            result.ref_end = last_anchor_end_ref + post_align.ref_end;
            temp_cigar_elements.insert(temp_cigar_elements.end(), post_align.cigar.begin(), post_align.cigar.end());
        }
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


const std::string GREEN = "\033[32m";
const std::string RED = "\033[31m";
const std::string YELLOW = "\033[33m";
const std::string BLUE = "\033[34m";
const std::string RESET = "\033[0m";

bool validate_test(const std::string& query, const std::string& reference, const std::vector<Anchor>& anchors, int k) {
    if (anchors.empty()) {
        std::cout << RED << "ERROR: No anchors provided" << RESET << std::endl;
        return false;
    }
    
    for (size_t i = 0; i < anchors.size(); ++i) {
        const Anchor& anchor = anchors[i];
        
        if (anchor.query_start + k > query.length()) {
            std::cout << RED << "ERROR: Anchor " << i << " query position " << anchor.query_start 
                      << " + k=" << k << " exceeds query length " << query.length() << RESET << std::endl;
            return false;
        }
        
        if (anchor.ref_start + k > reference.length()) {
            std::cout << RED << "ERROR: Anchor " << i << " ref position " << anchor.ref_start 
                      << " + k=" << k << " exceeds reference length " << reference.length() << RESET << std::endl;
            return false;
        }
        
        std::string query_kmer = query.substr(anchor.query_start, k);
        std::string ref_kmer = reference.substr(anchor.ref_start, k);
        
        if (query_kmer != ref_kmer) {
            std::cout << RED << "ERROR: Anchor " << i << " mismatch - Query: '" << query_kmer 
                      << "' vs Ref: '" << ref_kmer << "'" << RESET << std::endl;
            return false;
        }
    }
    
    for (size_t i = 1; i < anchors.size(); ++i) {
        if (anchors[i].query_start < anchors[i-1].query_start) {
            std::cout << RED << "ERROR: Anchors not in ascending query order" << RESET << std::endl;
            return false;
        }
        if (anchors[i].ref_start < anchors[i-1].ref_start) {
            std::cout << RED << "ERROR: Anchors not in ascending reference order" << RESET << std::endl;
            return false;
        }
    }
    
    return true;
}

bool validate_alignment(const std::string& query, const std::string& reference, const AlignmentResult& result) {
    if (result.query_start < 0 || result.query_end > query.length()) {
        std::cout << RED << "ERROR: Query alignment bounds invalid" << RESET << std::endl;
        return false;
    }
    
    if (result.ref_start < 0 || result.ref_end > reference.length()) {
        std::cout << RED << "ERROR: Reference alignment bounds invalid" << RESET << std::endl;
        return false;
    }
    
    if (result.query_start > result.query_end || result.ref_start > result.ref_end) {
        std::cout << RED << "ERROR: Alignment start > end" << RESET << std::endl;
        return false;
    }
    
    size_t query_pos = result.query_start;
    size_t ref_pos = result.ref_start;
    
    for (const auto& op : result.cigar) {
        switch (op.op) {
            case Operation::M:
            case Operation::Eq:
            case Operation::X:
                query_pos += op.len;
                ref_pos += op.len;
                break;
            case Operation::I:
                query_pos += op.len;
                break;
            case Operation::D:
                ref_pos += op.len;
                break;
            case Operation::Sentinel:
                break;
        }
        
        if (query_pos > query.length() || ref_pos > reference.length()) {
            std::cout << RED << "ERROR: CIGAR operations exceed sequence bounds" << RESET << std::endl;
            return false;
        }
    }
    
    if (query_pos != result.query_end || ref_pos != result.ref_end) {
        std::cout << RED << "ERROR: CIGAR end positions don't match alignment end positions" << RESET << std::endl;
        std::cout << "Expected query end: " << result.query_end << ", CIGAR end: " << query_pos << std::endl;
        std::cout << "Expected ref end: " << result.ref_end << ", CIGAR end: " << ref_pos << std::endl;
        return false;
    }
    
    return true;
}

struct TestCase {
    std::string name;
    std::string query;
    std::string reference;
    std::vector<Anchor> anchors;
    int k;
    int padding;
};

int main() {
    AlignmentScoring default_scoring = {
        .match = 3,
        .mismatch = -1,
        .gap_open = -3,
        .gap_extend = -1
    };

    std::vector<TestCase> test_cases = {
        // Test 1: Basic case (original)
        {
            "Basic alignment",
            "TCACTAACCGCTACGAT",
            "AAAATCACTACCCGCATACGTTCCCC",
            {{2, 6}, {7, 11}, {12, 17}},
            3, 2
        },
        
        // Test 2: Anchor at start of query
        {
            "Anchor at query start",
            "ATCGATCG",
            "AATCGGGGATCG",
            {{0, 1}, {5, 9}},
            3, 2
        },
        
        // Test 3: Anchor at start of reference
        {
            "Anchor at reference start",
            "GGGATCGATCG",
            "ATCGATCG",
            {{3, 0}, {6, 3}},
            3, 2
        },
        
        // Test 4: Anchor at end of query
        {
            "Anchor at query end",
            "GGGATCG",
            "AAGGATCGGGG",
            {{1, 2}, {4, 5}},
            3, 2
        },
        
        // Test 5: Anchor at end of reference
        {
            "Anchor at reference end",
            "AAAGATCGGGG",
            "GGGATCG",
            {{3, 2}, {5, 4}},
            3, 2
        },
        
        // Test 6: Single anchor
        {
            "Single anchor",
            "AAATCGAAA",
            "GGGGTCGGGGG",
            {{3, 4}},
            3, 2
        },
        
        // Test 7: Overlapping anchors (query)
        {
            "Overlapping anchors on query",
            "ATCGATCGATCG",
            "AAATCGAAAGATAAATCGAAA",
            {{1, 3}, {3, 9}, {5, 15}},
            3, 2
        },
        
        // Test 8: Overlapping anchors (reference)  
        {
            "Overlapping anchors on reference",
            "AAATCGAAAGATAAATCGAAA",
            "ATCGATCGATCG",
            {{3, 1}, {9, 3}, {16, 6}},
            3, 2
        },
        
        // Test 9: Adjacent anchors
        {
            "Adjacent anchors",
            "ATCGATCG",
            "ATCGATCG",
            {{0, 0}, {3, 3}},
            3, 2
        },
        
        // Test 10: Large gaps between anchors
        {
            "Large gaps between anchors",
            "ATCGAAAAAAAAAAGATCG",
            "ATCGGGGGGGGGGGGATCG",
            {{0, 0}, {16, 16}},
            3, 2
        },
        
        // Test 11: Different k value (k=4)
        {
            "Different k value (k=4)",
            "ATCGATCGATCG",
            "AAATCGAAATCGAAA",
            {{0, 2}, {4, 8}},
            4, 2
        },
        
        // Test 12: Different k value (k=2)
        {
            "Different k value (k=2)",
            "ATCGATCG",
            "AAATCGCGATCGAAA",
            {{0, 2}, {2, 6}, {4, 8}},
            2, 2
        },
        
        // Test 13: Many small anchors
        {
            "Many small anchors",
            "ATCGATCGATCG",
            "AATCGATCGATCGA",
            {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}},
            2, 1
        },
        
        // Test 14: Query much shorter than reference
        {
            "Short query, long reference",
            "ATCG",
            "GGGGGGGATCGGGGGGG",
            {{0, 7}},
            4, 3
        },
        
        // Test 15: Reference much shorter than query
        {
            "Long query, short reference",
            "GGGGGGGATCGGGGGGG",
            "ATCG",
            {{7, 0}},
            4, 3
        },
        
        // Test 16: Edge case - anchor at very end
        {
            "Anchor at sequence ends",
            "GGGATC",
            "GGGATC",
            {{3, 3}},
            3, 1
        },
        
        // Test 17: Complex overlapping pattern
        {
            "Complex overlapping pattern",
            "ATCGATCGATCGATCG",
            "GATCGATCGATCGATCGAA",
            {{0, 1}, {4, 5}, {8, 9}, {12, 13}},
            4, 2
        },
        
        // Test 18: Identical sequences
        {
            "Identical sequences",
            "ATCGATCG",
            "ATCGATCG",
            {{0, 0}, {4, 4}},
            4, 1
        },
        
        // Test 19: High padding value
        {
            "High padding value",
            "ATCGAAAAAGATCG",
            "ATCGGGGGGATCG",
            {{0, 0}, {10, 9}},
            4, 10
        },
        
        // Test 20: Multiple overlapping anchors
        {
            "Multiple overlapping anchors",
            "ATCGATCGATCGATCG",
            "AATCGATCGATCGATCGAA",
            {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, 13}},
            4, 2
        }
    };

    int passed_tests = 0;
    int total_tests = test_cases.size();

    std::cout << BLUE << "=== PIECEWISE ALIGNMENT TEST SUITE ===" << RESET << std::endl;
    std::cout << BLUE << "Running " << total_tests << " tests..." << RESET << std::endl << std::endl;

    for (size_t i = 0; i < test_cases.size(); ++i) {
        const auto& test = test_cases[i];
        
        std::cout << YELLOW << "Test " << (i + 1) << ": " << test.name << RESET << std::endl;
        std::cout << "Query: " << test.query << std::endl;
        std::cout << "Ref:   " << test.reference << std::endl;
        std::cout << "Anchors: ";
        for (const auto& anchor : test.anchors) {
            std::cout << "(" << anchor.query_start << "," << anchor.ref_start << ") ";
        }
        std::cout << "k=" << test.k << std::endl;

        bool test_valid = validate_test(test.query, test.reference, test.anchors, test.k);

        if (!test_valid) {
            std::cout << RED << "❌ TEST FAILED: Invalid input parameters" << RESET << std::endl;
            std::cout << std::endl;
            continue;
        }

        AlignmentResult result = piecewise_extension_alignment(
            test.query, test.reference, test.anchors, test.k, test.padding, default_scoring);

        std::cout << "Score: " << result.score << std::endl;
        std::cout << "Query Range: " << result.query_start << " - " << result.query_end << std::endl;
        std::cout << "Ref Range:   " << result.ref_start << " - " << result.ref_end << std::endl;
        std::cout << "CIGAR: " << result.to_cigar_string() << std::endl;

        bool alignment_valid = validate_alignment(test.query, test.reference, result);

        if (alignment_valid) {
            std::cout << GREEN << "✅ TEST PASSED" << RESET << std::endl;
            passed_tests++;
        } else {
            std::cout << RED << "❌ TEST FAILED: Invalid alignment result" << RESET << std::endl;
        }

        std::cout << std::endl;
        visualize_alignment(test.query, test.reference, result);
        std::cout << "----------------------------------------" << std::endl << std::endl;
    }

    std::cout << BLUE << "=== FINAL REPORT ===" << RESET << std::endl;
    std::cout << "Tests passed: " << GREEN << passed_tests << RESET << "/" << total_tests << std::endl;
    std::cout << "Tests failed: " << RED << (total_tests - passed_tests) << RESET << "/" << total_tests << std::endl;
    
    if (passed_tests == total_tests) {
        std::cout << GREEN << "ALL TESTS PASSED!" << RESET << std::endl;
    } else {
        std::cout << RED << "SOME TESTS FAILED" << RESET << std::endl;
    }

    return 0;
}
