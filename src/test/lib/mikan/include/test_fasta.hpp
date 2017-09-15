#ifndef MIKAN_TEST_FASTA_HPP_
#define MIKAN_TEST_FASTA_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mk_typedef.hpp"
#include "mk_input.hpp"

class TestFasta : public ::testing::Test {
protected:
    char *IFNAME1;
    char *IFNAME2;

    virtual void SetUp() {
        dfile = STRINGIZE(TEST_DATA_PATH);

        ifile1 = dfile;
        ifile2 = dfile;
        ifile1 += IFNAME1;
        ifile2 += IFNAME2;

        seqan::clear(mirna_ids);
        seqan::clear(mrna_ids);
        seqan::clear(mirna_seqs);
        seqan::clear(mrna_seqs);
    }

    void read_files() {
        coreInput.set_file_names(ifile1, ifile2);
        (void) coreInput.load_seq_from_file();
    }

    mikan::MKInput coreInput;

    mikan::TCharStr dfile;
    mikan::TCharStr ifile1;
    mikan::TCharStr ifile2;
    mikan::TCharSet mirna_ids;
    mikan::TCharSet mrna_ids;
    mikan::TRNASet mirna_seqs;
    mikan::TRNASet mrna_seqs;
};

#endif //MIKAN_TEST_FASTA_HPP_
