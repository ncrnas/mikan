#ifndef MIKANTEST_UTILS_HPP_
#define MIKANTEST_UTILS_HPP_

#include "seqan/seq_io.h"

int gtest_compare_two_files(seqan::CharString const &f1, seqan::CharString const &f2);

int gtest_compare_two_files2(seqan::CharString const &f1, seqan::CharString const &f2,
                             unsigned score_fld, unsigned round_dec, int ubound);

void comp_scores(std::string &score1, std::string &score2, unsigned round_dec, int ubound);

int gtest_read_file(seqan::CharString const &fname, std::vector<std::string> &lines);

void comp_two_rnas(seqan::RnaString const &rnaseq1, seqan::RnaString const &rnaseq2);

void split_line(std::string &line, std::vector<std::string> &flds);

#endif /* MIKANTEST_UTILS_HPP_ */