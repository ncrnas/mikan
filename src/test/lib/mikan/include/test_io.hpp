#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"

class TestIO : public ::testing::Test
{
protected:
    char *IFNAME1;
    char *IFNAME2;
    char *OFNAME1;
    char *OFNAME2;
    char *OMPATH;

    virtual void SetUp() {
        dfile = STRINGIZE(TEST_DATA_PATH);

        ifile1 = dfile;
        ifile2 = dfile;
        ifile1 += IFNAME1;
        ifile2 += IFNAME2;

        ompath = dfile;
        ompath += OMPATH;
        ofile1 = ompath;
        ofile2 = ompath;
        ofile1 += OFNAME1;
        ofile2 += OFNAME2;

        seqan::clear(mirna_ids);
        seqan::clear(mrna_ids);
        seqan::clear(mirna_seqs);
        seqan::clear(mrna_seqs);
    }

    virtual void TearDown() {
    }

    int fread_res;
    seqan::CharString dfile;
    seqan::CharString ifile1;
    seqan::CharString ifile2;
    seqan::CharString ompath;
    seqan::CharString ofile1;
    seqan::CharString ofile2;
    seqan::StringSet<seqan::CharString> mirna_ids;
    seqan::StringSet<seqan::CharString> mrna_ids;
    seqan::StringSet<seqan::RnaString> mirna_seqs;
    seqan::StringSet<seqan::RnaString> mrna_seqs;
};
