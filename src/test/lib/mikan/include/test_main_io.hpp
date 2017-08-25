#ifndef MIKAN_TEST_MAIN_IO_HPP_
#define MIKAN_TEST_MAIN_IO_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mk_typedef.hpp"
#include "mk_input.hpp"

class TestIOBase : public ::testing::Test {
protected:
    char *IFNAME1;
    char *IFNAME2;
    char *O1FNAME1;
    char *O1FNAME2;
    char *O2FNAME1;
    char *O2FNAME2;
    char *OMPATH;

    TestIOBase() {
        IFNAME1 = (char *) "";
        IFNAME2 = (char *) "";
        O1FNAME1 = (char *) "";
        O1FNAME2 = (char *) "";
        O2FNAME1 = (char *) "";
        O2FNAME2 = (char *) "";
        OMPATH = (char *) "";
    }

    virtual void SetUp() {
        seqan::CharString dfile = STRINGIZE(TEST_DATA_PATH);

        ifile1 = dfile;
        ifile2 = dfile;
        ifile1 += IFNAME1;
        ifile2 += IFNAME2;

        seqan::CharString ompath = dfile;
        ompath += OMPATH;
        o1file1 = ompath;
        o1file2 = ompath;
        o1file1 += O1FNAME1;
        o1file2 += O1FNAME2;

        o2file1 = ompath;
        o2file2 = ompath;
        o2file1 += O2FNAME1;
        o2file2 += O2FNAME2;

        argc = 5;
        argv[0] = (char *) "program";
        argv[1] = seqan::toCString(ifile1);
        argv[2] = seqan::toCString(ifile2);
        argv[3] = seqan::toCString(o2file1);
        argv[4] = seqan::toCString(o2file2);

        seqan::clear(mirna_ids);
        seqan::clear(mrna_ids);
        seqan::clear(mirna_seqs);
        seqan::clear(mrna_seqs);
    }

    void read_and_set_seqs() {
        read_files();
        set_seqs();
    }

    void read_files() {
        seqInput.set_file_names(ifile1, ifile2);
        (void) seqInput.load_seq_from_file();
    }

    void set_seqs() {
        mirna_ids = seqInput.get_mirna_ids();
        mirna_seqs = seqInput.get_mirna_seqs();
        mrna_ids = seqInput.get_mrna_ids();
        mrna_seqs = seqInput.get_mrna_seqs();
    }

    int argc;
    char *argv[6];

    mikan::MKInput seqInput;

    seqan::CharString ifile1;
    seqan::CharString ifile2;
    seqan::CharString o1file1;
    seqan::CharString o1file2;
    seqan::CharString o2file1;
    seqan::CharString o2file2;

    seqan::StringSet<seqan::CharString> mirna_ids;
    seqan::StringSet<seqan::CharString> mrna_ids;
    seqan::StringSet<seqan::RnaString> mirna_seqs;
    seqan::StringSet<seqan::RnaString> mrna_seqs;
};

#endif //MIKAN_TEST_MAIN_IO_HPP_
