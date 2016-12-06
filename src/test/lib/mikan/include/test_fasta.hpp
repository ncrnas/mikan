#ifndef MIKAN_TEST_FASTA_HPP_
#define MIKAN_TEST_FASTA_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mr3_core.hpp"

class TestFasta : public ::testing::Test
{
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
        options.mMiRNAFasta = seqan::toCString(ifile1);
        options.mMRNAFasta = seqan::toCString(ifile2);
        coreInput.init_from_args(options);
        (void)coreInput.load_seq_from_file();
    }

    mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
    mr3as::MR3Options options;

    seqan::CharString dfile;
    seqan::CharString ifile1;
    seqan::CharString ifile2;
    seqan::StringSet<seqan::CharString> mirna_ids;
    seqan::StringSet<seqan::CharString> mrna_ids;
    seqan::StringSet<seqan::RnaString> mirna_seqs;
    seqan::StringSet<seqan::RnaString> mrna_seqs;
};

#endif //MIKAN_TEST_FASTA_HPP_
