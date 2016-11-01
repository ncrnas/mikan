#include "gtest/gtest.h"
#include <seqan/seq_io.h>
#include "mikan_utils.hpp"

int gtest_compare_two_files(seqan::CharString const &f1, seqan::CharString const &f2) {
    seqan::StringSet<seqan::CharString> lines1;
    seqan::StringSet<seqan::CharString> lines2;
    int fres1, fres2;

    fres1 = gtest_read_file(f1, lines1);
    if (fres1 != 0 )
    {
        std::cerr << "ERROR: Could not open file1" << std::endl;
        return 1;
    }

    fres2 = gtest_read_file(f2, lines2);
    if (fres2 != 0)
    {
        std::cerr << "ERROR: Could not open file2" << std::endl;
        return 1;
    }

    EXPECT_EQ(length(lines1), length(lines2));

    for (unsigned i = 0; i < length(lines1); ++i) {
        EXPECT_STREQ(seqan::toCString((seqan::CharString)lines1[i]),
                     seqan::toCString((seqan::CharString)lines2[i]));
    }

    return 0;
}

int gtest_read_file(seqan::CharString const &fname,
                    seqan::StringSet<seqan::CharString> &lines)
{
    seqan::CharString line;

    std::ifstream file(seqan::toCString(fname));
    std::string str;
    while (std::getline(file, str)) {
        line = str;
        appendValue(lines, line);
    }

    return 0;
}

void comp_two_rnas(seqan::RnaString const &rnaseq1, seqan::RnaString const &rnaseq2)
{
    EXPECT_STREQ(seqan::toCString((seqan::CharString)rnaseq1),
                 seqan::toCString((seqan::CharString)rnaseq2));
}
