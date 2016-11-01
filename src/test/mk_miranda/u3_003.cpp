#include <iostream>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class U3003 : public TestIO
    {
    protected:
        U3003() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_003.fasta";
            OFNAME1 = (char *)"test_output_1.txt";
            OFNAME2 = (char *)"test_ds1.txt";
            OMPATH = (char *)"miranda/";
        }

        void read_files() {
            options.mMiRNAFasta = ifile1;
            options.mMRNAFasta = ifile2;

            coreInput.init_from_args(options);
            int retVal = coreInput.load_seq_from_file();
            EXPECT_EQ(retVal, 0);
        }
        mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
        mr3as::MR3Options options;
    };

    TEST_F(U3003, mrna_fasta) {
        read_files();
        mrna_ids = coreInput.get_mrna_ids();
        mrna_seqs = coreInput.get_mrna_seqs();

        EXPECT_EQ(1u, length(mrna_ids));

        const char *id1 = "hg18_refGene NM_018214 range=chr6:53895551-53896878 5'pad=0 3'pad=0 "
                "revComp=FALSE strand=+ repeatMasking=none";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *seq1 = "AGUUUCACCUCCAAGUUUUACCUCCUGUGUCUUCCUCUGCUGUCGAGACG"
                "UUCCUGUCUGCUUCCCGGGAGCCUCACGUGCUCCUUGUCCUAACCAGCCC"
                "CCGCGCGCCAUCUUCCCGUGGAGUGUGGGGAAGCUGCUGUCUCCCAGGAA"
                "GUGCCUUACUCAUCCCGCAACCAGUCAGCGCACCAGUGGUCUCCCGGUGU"
                "GAUUUUUUUUUUUUUUAAUUUCAGUUGUUUGUAAUAAGUAGAAUACACUA"
                "CUGUAAACAUACGACCUUUGUUUUUGUCUUAUGUUGGGGUAAAGGAAAGC"
                "AGGAAGGGGAAUUUUUAUCCUCCUCCCUUCCGUAAAGUGCUGGGAUAUUU"
                "UGAAUCCCCCAAGUUCCCUUGGACCUACUGAUGAGAGAUAGUUUUAUGUA"
                "UGGGGAAAAAUGGAUACUUUUUAAACCUUUUUUGGCAGCUCAGAUGGUGU"
                "AAAUUUUAAAAUUUUGUAUAGGUAUUUCAUAACAAAAAUAUGUAUUUCUU"
                "UUUUGUUAUUUUAUCUUGAAAACGGUACAUAUUUUAGUAUUUGUGCAGAA"
                "AAACAAGUCCUAAAGUAUUUGUUUUUAUUUGUACCAUCCACUUGUGCCUU"
                "ACUGUAUCCUGUGUCAUGUCCAAUCAGUUGUAAACAAUGGCAUCUUUGAA"
                "CAGUGUGAUGAGAAUAGGAAUGUGGUGUUUUAAAGCAGUGUUGCAUUUUA"
                "AUCAGUAAUCUACCUGGUGGAUUUGUUUUUAACCAAAAAGAUGAAUUAUC"
                "AAUGAUUUGUAAUUAUAUCAGUUGAUUUUUUUUGAAAAGAUGAACCAAAG"
                "GAUUUGACUGCUAAUAUUUUAUUCCUUACACUUUUUUUCUGAAUAAGUCU"
                "CUCAUAAUGAGUGCAGUGUCAGACUGUGCCUACUCUGAUGGUAUGUGCCA"
                "UUUGUAAAAUAAAAUAGAGCAGAAAAACACAAAAAGAGAACACUGGUUCA"
                "GACAUUCAGUGGGCAAGUAAAUUAUGGACUGCAAAAUAAUGAUUUUUAUU"
                "CAAGAAAGCUUUAAAAGUUUUAUAUCCAGAUAUACAACCACAAUAAAGCA"
                "AAAUAACCUACUAUCAAAAUAGAAAUGUUGCUAUCUUUAUAAGUGCAAUU"
                "UAAUUUGUAAAUAGAGUUUGAAUCAAAGUAUCACAAAAUACUGCUUCAAG"
                "AUUUAAUUUUAAAUCUGCUAAUUUAAGGGAUAUUGGGAAAAGUUUUGGUG"
                "UGUUUCUGUUGAUUUCUUUUUUGUAUGCUGUGAUAAAAGAGAAAUGAAAA"
                "GUGCCAGUCACUGUGUGGUGUCUAGGAAAAUCAUAUAUAUUUUUUUCUCC"
                "AAGAAAUAAAUUCAUCCUGGACAUUGGC";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }
}