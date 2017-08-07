#ifndef MIKANTEST_UTILS_HPP_
#define MIKANTEST_UTILS_HPP_

#include "seqan/seq_io.h"

int gtest_compare_two_files(seqan::CharString const &f1, seqan::CharString const &f2);

int gtest_compare_two_files2(seqan::CharString const &f1, seqan::CharString const &f2,
                             unsigned score_fld, unsigned round_dec, int ubound);

int gtest_compare_two_files3(seqan::CharString const &f1, seqan::CharString const &f2,
                             bool uppercase, bool replace_space);

void comp_two_rnas(seqan::RnaString const &rnaseq1, seqan::RnaString const &rnaseq2);


#endif /* MIKANTEST_UTILS_HPP_ */