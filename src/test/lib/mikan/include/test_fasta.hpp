#ifndef MIKAN_TEST_FASTA_HPP_
#define MIKAN_TEST_FASTA_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mk_inst_template.hpp"
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

    mikan::MKInput<mikan::TRNATYPE> coreInput;
    
    seqan::CharString dfile;
    seqan::CharString ifile1;
    seqan::CharString ifile2;
    seqan::StringSet<seqan::CharString> mirna_ids;
    seqan::StringSet<seqan::CharString> mrna_ids;
    seqan::StringSet<seqan::RnaString> mirna_seqs;
    seqan::StringSet<seqan::RnaString> mrna_seqs;
};

#endif //MIKAN_TEST_FASTA_HPP_
