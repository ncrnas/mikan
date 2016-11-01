#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mr3_core.hpp"

class TestIOCommon : public ::testing::Test
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

class TestIOBase : public ::testing::Test
{
protected:
    char *IFNAME1;
    char *IFNAME2;
    char *O1FNAME1;
    char *O1FNAME2;
    char *O2FNAME1;
    char *O2FNAME2;
    char *OMPATH;

    virtual void SetUp() {
        dfile = STRINGIZE(TEST_DATA_PATH);

        ifile1 = dfile;
        ifile2 = dfile;
        ifile1 += IFNAME1;
        ifile2 += IFNAME2;

        ompath = dfile;
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
        argv[0] = (char *)"program";
        argv[1] = seqan::toCString(ifile1);
        argv[2] = seqan::toCString(ifile2);
        argv[3] = seqan::toCString(o2file1);
        argv[4] = seqan::toCString(o2file2);

        seqan::clear(mirna_ids);
        seqan::clear(mrna_ids);
        seqan::clear(mirna_seqs);
        seqan::clear(mrna_seqs);
    }

    int fread_res;
    int argc;
    char *argv[5];
    seqan::CharString dfile;
    seqan::CharString ifile1;
    seqan::CharString ifile2;
    seqan::CharString ompath;
    seqan::CharString o1file1;
    seqan::CharString o1file2;
    seqan::CharString o2file1;
    seqan::CharString o2file2;
    seqan::StringSet<seqan::CharString> mirna_ids;
    seqan::StringSet<seqan::CharString> mrna_ids;
    seqan::StringSet<seqan::RnaString> mirna_seqs;
    seqan::StringSet<seqan::RnaString> mrna_seqs;
};

class TestIOMR3AS : public TestIOBase
{
protected:

    void read_files(bool parse_argv) {
        if (parse_argv)
        {
            (void)options.parseCommandLine(argc, (const char **)argv);
        }
        else
        {
            options.mMiRNAFasta = seqan::toCString(ifile1);
            options.mMRNAFasta = seqan::toCString(ifile2);
            coreInput.init_from_args(options);
            (void)coreInput.load_seq_from_file();
        }

        coreInput.init_from_args(options);
        (void)coreInput.load_seq_from_file();
    }

    void run_main() {
        (void)mr3as::MR3CoreMain(argc, (const char **)argv);
    }

    mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
    mr3as::MR3Options options;
    mr3as::MR3SeedSeqs<seqan::RnaString> mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedDef;
};