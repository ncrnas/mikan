#include <string.h>
#include <algorithm>
#include "gtest/gtest.h"
#include <seqan/seq_io.h>
#include "mikan_utils.hpp"

// Static functions
static int gtest_read_files(seqan::CharString const &fname1, seqan::CharString const &fname2,
                            std::vector<std::string> &lines1, std::vector<std::string> &lines2);

static int gtest_read_file(seqan::CharString const &fname, std::vector<std::string> &lines);

static void comp_scores(std::string &score1, std::string &score2, unsigned round_dec, int ubound);

static void split_line(std::string &line, std::vector<std::string> &flds);

int gtest_compare_two_files(seqan::CharString const &f1, seqan::CharString const &f2) {
    std::vector<std::string> lines1;
    std::vector<std::string> lines2;

    int res = gtest_read_files(f1, f2, lines1, lines2);
    if (res != 0) {
        return 1;
    }

    for (unsigned i = 0; i < lines1.size(); ++i) {
        EXPECT_STREQ(lines1[i].c_str(), lines2[i].c_str());
    }

    return 0;
}

int gtest_compare_two_files2(
        seqan::CharString const &f1,
        seqan::CharString const &f2,
        unsigned score_fld,
        unsigned round_dec,
        int ubound) {

    std::vector<std::string> lines1;
    std::vector<std::string> lines2;

    int res = gtest_read_files(f1, f2, lines1, lines2);
    if (res != 0) {
        return 1;
    }

    for (unsigned i = 0; i < lines1.size(); ++i) {
        std::vector<std::string> flds1;
        std::vector<std::string> flds2;
        split_line(lines1[i], flds1);
        split_line(lines2[i], flds2);

//        std::cout << lines1[i] << std::endl;

        EXPECT_EQ(flds1.size(), flds2.size());

        for (unsigned j = 0; j < flds1.size(); ++j) {
            if (score_fld != 0 && j == score_fld) {
                comp_scores(flds1[j], flds2[j], round_dec, ubound);
            } else {
//                std::cout << flds1[j] << " == " << flds2[j] << std::endl;
                EXPECT_STREQ(flds1[j].c_str(), flds2[j].c_str());
            }
        }

    }

    return 0;
}

int gtest_compare_two_files3(
        seqan::CharString const &f1,
        seqan::CharString const &f2,
        bool uppercase,
        bool replace_space) {

    std::vector<std::string> lines1;
    std::vector<std::string> lines2;

    int res = gtest_read_files(f1, f2, lines1, lines2);
    if (res != 0) {
        return 1;
    }

    for (unsigned i = 0; i < lines1.size(); ++i) {
        std::vector<std::string> flds1;
        std::vector<std::string> flds2;
        split_line(lines1[i], flds1);
        split_line(lines2[i], flds2);

        for (unsigned j = 0; j < 2; ++j) {
            if (uppercase) {
                std::transform(flds1[j].begin(), flds1[j].end(), flds1[j].begin(), toupper);
                std::transform(flds2[j].begin(), flds2[j].end(), flds2[j].begin(), toupper);
            }

            if (replace_space) {
                std::replace(flds1[j].begin(), flds1[j].end(), ' ', '_');
                std::replace(flds2[j].begin(), flds2[j].end(), ' ', '_');
            }

            EXPECT_STREQ(flds1[j].c_str(), flds2[j].c_str());
        }

    }

    return 0;
}

void comp_two_rnas(seqan::RnaString const &rnaseq1, seqan::RnaString const &rnaseq2) {
    EXPECT_STREQ(seqan::toCString((seqan::CharString) rnaseq1),
                 seqan::toCString((seqan::CharString) rnaseq2));
}

static int gtest_read_files(seqan::CharString const &fname1,
                     seqan::CharString const &fname2,
                     std::vector<std::string> &lines1,
                     std::vector<std::string> &lines2) {

    int fres1, fres2;

    fres1 = gtest_read_file(fname1, lines1);
    if (fres1 != 0) {
        std::cerr << "ERROR: Could not open file1" << std::endl;
        return 1;
    }

    fres2 = gtest_read_file(fname2, lines2);
    if (fres2 != 0) {
        std::cerr << "ERROR: Could not open file2" << std::endl;
        return 1;
    }

    std::sort(lines1.begin(), lines1.end());
    std::sort(lines2.begin(), lines2.end());

    EXPECT_EQ(lines1.size(), lines2.size());

    return 0;
}


static int gtest_read_file(seqan::CharString const &fname,
                    std::vector<std::string> &lines) {

    std::ifstream file(seqan::toCString(fname));
    std::string str;

    while (std::getline(file, str)) {
        lines.push_back(str);
    }

    return 0;
}

static void split_line(std::string &line, std::vector<std::string> &flds) {
    char *pline = const_cast<char *>(line.c_str());
    char *pch;

    flds.clear();

    pch = strtok(pline, "\t");

    while (pch != NULL) {
        std::string str(pch);
        if (str.size() > 0) {
            flds.push_back(str);
        }
        pch = strtok(NULL, "\t");
    }

}

static void comp_scores(std::string &score1, std::string &score2, unsigned round_dec, int ubound) {

    double ds1 = ::atof(score1.c_str());
    double ds2 = ::atof(score2.c_str());
    int s1 = static_cast<int>(ds1 * round_dec);
    int s2 = static_cast<int>(ds2 * round_dec);

    int diff  = s1 - s2;
    int lbound = -1 * ubound;

//    std::cout << s1 << " == " << s2  << std::endl;
//    std::cout << lbound << " <= " << diff << " <= " << ubound << std::endl;
//    EXPECT_EQ(s1, s2);

    EXPECT_TRUE(lbound <= diff && diff <= ubound);

}
