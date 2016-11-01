#include <iostream>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class U3005 : public TestIOCommon
    {
    protected:
        U3005() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_005.fasta";
        }
    };

    TEST_F(U3005, mrna_fasta) {
        read_files();
        mrna_ids = coreInput.get_mrna_ids();
        mrna_seqs = coreInput.get_mrna_seqs();

        EXPECT_EQ(6u, length(mrna_ids));

        const char *id1 = "hg18_refgene test1_8mer";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        const char *id2 = "hg18_refgene test2_7mer-m8";
        EXPECT_STREQ(id2, seqan::toCString(mrna_ids[1]));

        const char *id3 = "hg18_refgene test3_7mer-a1";
        EXPECT_STREQ(id3, seqan::toCString(mrna_ids[2]));

        const char *id4 = "hg18_refgene test4_6mer1";
        EXPECT_STREQ(id4, seqan::toCString(mrna_ids[3]));

        const char *id5 = "hg18_refgene test5_6mer2";
        EXPECT_STREQ(id5, seqan::toCString(mrna_ids[4]));

        const char *id6 = "hg18_refgene test5_6mer3";
        EXPECT_STREQ(id6, seqan::toCString(mrna_ids[5]));

        EXPECT_EQ(6u, length(mrna_seqs));
        const char *seq1 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCCCAGUGCCUUACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCUGGAAAAAACCUUAAAAAAUGUUUCCUCCAAAUCUGAUUUCAU"
                "UACAUUUCUGAAUUGUUGGGGUUUUUUUUGUUGUUUUGUUUUGUUUUGUA"
                "GAUGGAGUUUCACUUUUGUUGCCCAGGCUGGAGUGUAGUGGCGCGAUCUC"
                "GGCUCAGCCUCCCGAGUAGCUGGGAUUACAGGCAUGUGCCACCACGCCCG"
                "GCUAAUGUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUGGUCAGG"
                "CUGGUCUCAAACUCCUGACCUCAGGUGAUCCGCCCACCUCAGCCUCCCAA"
                "AGUGCUGGGAUGACAGGUGUGAGCCACUGCGCCCAGCCUGAAUCAUUUCU"
                "UAUACCUUCUGACAGCCCAACUUCCAGAGGACAGCUCUGGGGUACUCGUU"
                "GGAUGUCUGUGAGUACCUGGUCAUACGGGUCAGUAGGGAUAAGAAUUGUC"
                "UCUGGGCUGAGGAAUUCUUCUGUUCUCUGGUUUCACCAGCGUUGGGUUUG"
                "CUCAUGUAAUGUGGUCACCAUACUCAAAUGGUGUCAUGGCUGAAGUUGGC"
                "CACCUUGCUUGAGGGACAAGUUGUUUAUGUAUCAGCUCUCUGCUGGGUCU"
                "CCCUUUCCAUGGCAAAUGGGCAGCUCCAUCCUCUUGACUCUUCUAAAUGC"
                "CCAAAAGAGGUGUCAUGCUUUGGGGGUACGAUGUUUAUACUCCGUAAAGA"
                "ACAUACAAGGACAUUCACUGCUGAUUUUUUUUUUUGUUUGUUUGAGACAG"
                "GGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGAUGCAAUCUUGGCUCA"
                "CUGCAACCUCCGCCUCUCAGGUUCAAGUGGUUCUCCUUCCUCAGCCUCCC"
                "AAGUAGCUGGGAUUACAGGCACCUACCACCAGGGCCAGCUAAUUUUUGUA"
                "UGUUUAGUAGAAACGGGGUUUCACCAUGUUGGCCAGGCUGUUCUCGAACU"
                "CCUGACCUCAGGUGAUCUGCCCGCCUCGGUCUCCCAAAGUGCUGGGAUUA"
                "CAGGCAUGAGCCACUGCACCUGACCUGCUGAAUUGUUUAUAAUGGCAAGA"
                "AAUAGGAAACCCCCCAAUGUCUGUUGAACAGCUAUCACGUUGAACCACGU"
                "GAAACUGCUGUUUUCUAGGCCAAAAAUGGUGAGCGAUCAUUUAUUUCAUG"
                "AUUCAACCUGAUACAUUUACAUAGUGCAAAACUGUGUCACAGUUUCAGGC"
                "UUUUAUGAGGAAAGCGUUUCUGUGUAGAAACUGGAAGCUGUUCAGGGCAU"
                "CGGCAGCUGAACCCUGCUCCGUUGGUCAGCGUUACUAUCAUCUCGGAUCA"
                "UAUGGAGCUCAUGUCAGCCGUGUGGGUGGCGGGUGCACAGAGACGGUCUG"
                "GAAGGAAACACGCGGAUCUGAACAGCAGUAAUCCUGGGGGAUACGGGGGU"
                "UGGGCUAGAUUACAGAGGGCUCAUUUUCUACGUCAUGUAUUUUAUGAUAC"
                "UUGAAUUUUUUGAAAUGGGCAUUUAUUUUAUAACAUGUUAAAAUGUACUU"
                "UUUAAAUUAAGUCAUUUUGUAAUAUUUGAAUUUUUACAUUUGUUGUACAA"
                "UCAGGAAAAGCAAUAAAGAUUUUUCAAAAAUAG";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));

        const char *seq2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCCCAGUGUCUUCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCUGGAAAAAACCUUAAAAAAUGUUUCCUCCAAAUCUGAUUUCAU"
                "UACAUUUCUGAAUUGUUGGGGUUUUUUUUGUUGUUUUGUUUUGUUUUGUA"
                "GAUGGAGUUUCACUUUUGUUGCCCAGGCUGGAGUGUAGUGGCGCGAUCUC"
                "GGCUCAGCCUCCCGAGUAGCUGGGAUUACAGGCAUGUGCCACCACGCCCG"
                "GCUAAUGUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUGGUCAGG"
                "CUGGUCUCAAACUCCUGACCUCAGGUGAUCCGCCCACCUCAGCCUCCCAA"
                "AGUGCUGGGAUGACAGGUGUGAGCCACUGCGCCCAGCCUGAAUCAUUUCU"
                "UAUACCUUCUGACAGCCCAACUUCCAGAGGACAGCUCUGGGGUACUCGUU"
                "GGAUGUCUGUGAGUACCUGGUCAUACGGGUCAGUAGGGAUAAGAAUUGUC"
                "UCUGGGCUGAGGAAUUCUUCUGUUCUCUGGUUUCACCAGCGUUGGGUUUG"
                "CUCAUGUAAUGUGGUCACCAUACUCAAAUGGUGUCAUGGCUGAAGUUGGC"
                "CACCUUGCUUGAGGGACAAGUUGUUUAUGUAUCAGCUCUCUGCUGGGUCU"
                "CCCUUUCCAUGGCAAAUGGGCAGCUCCAUCCUCUUGACUCUUCUAAAUGC"
                "CCAAAAGAGGUGUCAUGCUUUGGGGGUACGAUGUUUAUACUCCGUAAAGA"
                "ACAUACAAGGACAUUCACUGCUGAUUUUUUUUUUUGUUUGUUUGAGACAG"
                "GGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGAUGCAAUCUUGGCUCA"
                "CUGCAACCUCCGCCUCUCAGGUUCAAGUGGUUCUCCUUCCUCAGCCUCCC"
                "AAGUAGCUGGGAUUACAGGCACCUACCACCAGGGCCAGCUAAUUUUUGUA"
                "UGUUUAGUAGAAACGGGGUUUCACCAUGUUGGCCAGGCUGUUCUCGAACU"
                "CCUGACCUCAGGUGAUCUGCCCGCCUCGGUCUCCCAAAGUGCUGGGAUUA"
                "CAGGCAUGAGCCACUGCACCUGACCUGCUGAAUUGUUUAUAAUGGCAAGA"
                "AAUAGGAAACCCCCCAAUGUCUGUUGAACAGCUAUCACGUUGAACCACGU"
                "GAAACUGCUGUUUUCUAGGCCAAAAAUGGUGAGCGAUCAUUUAUUUCAUG"
                "AUUCAACCUGAUACAUUUACAUAGUGCAAAACUGUGUCACAGUUUCAGGC"
                "UUUUAUGAGGAAAGCGUUUCUGUGUAGAAACUGGAAGCUGUUCAGGGCAU"
                "CGGCAGCUGAACCCUGCUCCGUUGGUCAGCGUUACUAUCAUCUCGGAUCA"
                "UAUGGAGCUCAUGUCAGCCGUGUGGGUGGCGGGUGCACAGAGACGGUCUG"
                "GAAGGAAACACGCGGAUCUGAACAGCAGUAAUCCUGGGGGAUACGGGGGU"
                "UGGGCUAGAUUACAGAGGGCUCAUUUUCUACGUCAUGUAUUUUAUGAUAC"
                "UUGAAUUUUUUGAAAUGGGCAUUUAUUUUAUAACAUGUUAAAAUGUACUU"
                "UUUAAAUUAAGUCAUUUUGUAAUAUUUGAAUUUUUACAUUUGUUGUACAA"
                "UCAGGAAAAGCAAUAAAGAUUUUUCAAAAAUAG";
        EXPECT_STREQ(seq2, seqan::toCString((seqan::CharString)mrna_seqs[1]));

        const char *seq3 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCCCAAUGUCUUACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCUGGAAAAAACCUUAAAAAAUGUUUCCUCCAAAUCUGAUUUCAU"
                "UACAUUUCUGAAUUGUUGGGGUUUUUUUUGUUGUUUUGUUUUGUUUUGUA"
                "GAUGGAGUUUCACUUUUGUUGCCCAGGCUGGAGUGUAGUGGCGCGAUCUC"
                "GGCUCAGCCUCCCGAGUAGCUGGGAUUACAGGCAUGUGCCACCACGCCCG"
                "GCUAAUGUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUGGUCAGG"
                "CUGGUCUCAAACUCCUGACCUCAGGUGAUCCGCCCACCUCAGCCUCCCAA"
                "AGUGCUGGGAUGACAGGUGUGAGCCACUGCGCCCAGCCUGAAUCAUUUCU"
                "UAUACCUUCUGACAGCCCAACUUCCAGAGGACAGCUCUGGGGUACUCGUU"
                "GGAUGUCUGUGAGUACCUGGUCAUACGGGUCAGUAGGGAUAAGAAUUGUC"
                "UCUGGGCUGAGGAAUUCUUCUGUUCUCUGGUUUCACCAGCGUUGGGUUUG"
                "CUCAUGUAAUGUGGUCACCAUACUCAAAUGGUGUCAUGGCUGAAGUUGGC"
                "CACCUUGCUUGAGGGACAAGUUGUUUAUGUAUCAGCUCUCUGCUGGGUCU"
                "CCCUUUCCAUGGCAAAUGGGCAGCUCCAUCCUCUUGACUCUUCUAAAUGC"
                "CCAAAAGAGGUGUCAUGCUUUGGGGGUACGAUGUUUAUACUCCGUAAAGA"
                "ACAUACAAGGACAUUCACUGCUGAUUUUUUUUUUUGUUUGUUUGAGACAG"
                "GGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGAUGCAAUCUUGGCUCA"
                "CUGCAACCUCCGCCUCUCAGGUUCAAGUGGUUCUCCUUCCUCAGCCUCCC"
                "AAGUAGCUGGGAUUACAGGCACCUACCACCAGGGCCAGCUAAUUUUUGUA"
                "UGUUUAGUAGAAACGGGGUUUCACCAUGUUGGCCAGGCUGUUCUCGAACU"
                "CCUGACCUCAGGUGAUCUGCCCGCCUCGGUCUCCCAAAGUGCUGGGAUUA"
                "CAGGCAUGAGCCACUGCACCUGACCUGCUGAAUUGUUUAUAAUGGCAAGA"
                "AAUAGGAAACCCCCCAAUGUCUGUUGAACAGCUAUCACGUUGAACCACGU"
                "GAAACUGCUGUUUUCUAGGCCAAAAAUGGUGAGCGAUCAUUUAUUUCAUG"
                "AUUCAACCUGAUACAUUUACAUAGUGCAAAACUGUGUCACAGUUUCAGGC"
                "UUUUAUGAGGAAAGCGUUUCUGUGUAGAAACUGGAAGCUGUUCAGGGCAU"
                "CGGCAGCUGAACCCUGCUCCGUUGGUCAGCGUUACUAUCAUCUCGGAUCA"
                "UAUGGAGCUCAUGUCAGCCGUGUGGGUGGCGGGUGCACAGAGACGGUCUG"
                "GAAGGAAACACGCGGAUCUGAACAGCAGUAAUCCUGGGGGAUACGGGGGU"
                "UGGGCUAGAUUACAGAGGGCUCAUUUUCUACGUCAUGUAUUUUAUGAUAC"
                "UUGAAUUUUUUGAAAUGGGCAUUUAUUUUAUAACAUGUUAAAAUGUACUU"
                "UUUAAAUUAAGUCAUUUUGUAAUAUUUGAAUUUUUACAUUUGUUGUACAA"
                "UCAGGAAAAGCAAUAAAGAUUUUUCAAAAAUAG";
        EXPECT_STREQ(seq3, seqan::toCString((seqan::CharString)mrna_seqs[2]));

        const char *seq4 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCCCAGUGCCUCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCUGGAAAAAACCUUAAAAAAUGUUUCCUCCAAAUCUGAUUUCAU"
                "UACAUUUCUGAAUUGUUGGGGUUUUUUUUGUUGUUUUGUUUUGUUUUGUA"
                "GAUGGAGUUUCACUUUUGUUGCCCAGGCUGGAGUGUAGUGGCGCGAUCUC"
                "GGCUCAGCCUCCCGAGUAGCUGGGAUUACAGGCAUGUGCCACCACGCCCG"
                "GCUAAUGUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUGGUCAGG"
                "CUGGUCUCAAACUCCUGACCUCAGGUGAUCCGCCCACCUCAGCCUCCCAA"
                "AGUGCUGGGAUGACAGGUGUGAGCCACUGCGCCCAGCCUGAAUCAUUUCU"
                "UAUACCUUCUGACAGCCCAACUUCCAGAGGACAGCUCUGGGGUACUCGUU"
                "GGAUGUCUGUGAGUACCUGGUCAUACGGGUCAGUAGGGAUAAGAAUUGUC"
                "UCUGGGCUGAGGAAUUCUUCUGUUCUCUGGUUUCACCAGCGUUGGGUUUG"
                "CUCAUGUAAUGUGGUCACCAUACUCAAAUGGUGUCAUGGCUGAAGUUGGC"
                "CACCUUGCUUGAGGGACAAGUUGUUUAUGUAUCAGCUCUCUGCUGGGUCU"
                "CCCUUUCCAUGGCAAAUGGGCAGCUCCAUCCUCUUGACUCUUCUAAAUGC"
                "CCAAAAGAGGUGUCAUGCUUUGGGGGUACGAUGUUUAUACUCCGUAAAGA"
                "ACAUACAAGGACAUUCACUGCUGAUUUUUUUUUUUGUUUGUUUGAGACAG"
                "GGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGAUGCAAUCUUGGCUCA"
                "CUGCAACCUCCGCCUCUCAGGUUCAAGUGGUUCUCCUUCCUCAGCCUCCC"
                "AAGUAGCUGGGAUUACAGGCACCUACCACCAGGGCCAGCUAAUUUUUGUA"
                "UGUUUAGUAGAAACGGGGUUUCACCAUGUUGGCCAGGCUGUUCUCGAACU"
                "CCUGACCUCAGGUGAUCUGCCCGCCUCGGUCUCCCAAAGUGCUGGGAUUA"
                "CAGGCAUGAGCCACUGCACCUGACCUGCUGAAUUGUUUAUAAUGGCAAGA"
                "AAUAGGAAACCCCCCAAUGUCUGUUGAACAGCUAUCACGUUGAACCACGU"
                "GAAACUGCUGUUUUCUAGGCCAAAAAUGGUGAGCGAUCAUUUAUUUCAUG"
                "AUUCAACCUGAUACAUUUACAUAGUGCAAAACUGUGUCACAGUUUCAGGC"
                "UUUUAUGAGGAAAGCGUUUCUGUGUAGAAACUGGAAGCUGUUCAGGGCAU"
                "CGGCAGCUGAACCCUGCUCCGUUGGUCAGCGUUACUAUCAUCUCGGAUCA"
                "UAUGGAGCUCAUGUCAGCCGUGUGGGUGGCGGGUGCACAGAGACGGUCUG"
                "GAAGGAAACACGCGGAUCUGAACAGCAGUAAUCCUGGGGGAUACGGGGGU"
                "UGGGCUAGAUUACAGAGGGCUCAUUUUCUACGUCAUGUAUUUUAUGAUAC"
                "UUGAAUUUUUUGAAAUGGGCAUUUAUUUUAUAACAUGUUAAAAUGUACUU"
                "UUUAAAUUAAGUCAUUUUGUAAUAUUUGAAUUUUUACAUUUGUUGUACAA"
                "UCAGGAAAAGCAAUAAAGAUUUUUCAAAAAUAG";
        EXPECT_STREQ(seq4, seqan::toCString((seqan::CharString)mrna_seqs[3]));

        const char *seq5 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCCCAGUGUCUGUCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCUGGAAAAAACCUUAAAAAAUGUUUCCUCCAAAUCUGAUUUCAU"
                "UACAUUUCUGAAUUGUUGGGGUUUUUUUUGUUGUUUUGUUUUGUUUUGUA"
                "GAUGGAGUUUCACUUUUGUUGCCCAGGCUGGAGUGUAGUGGCGCGAUCUC"
                "GGCUCAGCCUCCCGAGUAGCUGGGAUUACAGGCAUGUGCCACCACGCCCG"
                "GCUAAUGUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUGGUCAGG"
                "CUGGUCUCAAACUCCUGACCUCAGGUGAUCCGCCCACCUCAGCCUCCCAA"
                "AGUGCUGGGAUGACAGGUGUGAGCCACUGCGCCCAGCCUGAAUCAUUUCU"
                "UAUACCUUCUGACAGCCCAACUUCCAGAGGACAGCUCUGGGGUACUCGUU"
                "GGAUGUCUGUGAGUACCUGGUCAUACGGGUCAGUAGGGAUAAGAAUUGUC"
                "UCUGGGCUGAGGAAUUCUUCUGUUCUCUGGUUUCACCAGCGUUGGGUUUG"
                "CUCAUGUAAUGUGGUCACCAUACUCAAAUGGUGUCAUGGCUGAAGUUGGC"
                "CACCUUGCUUGAGGGACAAGUUGUUUAUGUAUCAGCUCUCUGCUGGGUCU"
                "CCCUUUCCAUGGCAAAUGGGCAGCUCCAUCCUCUUGACUCUUCUAAAUGC"
                "CCAAAAGAGGUGUCAUGCUUUGGGGGUACGAUGUUUAUACUCCGUAAAGA"
                "ACAUACAAGGACAUUCACUGCUGAUUUUUUUUUUUGUUUGUUUGAGACAG"
                "GGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGAUGCAAUCUUGGCUCA"
                "CUGCAACCUCCGCCUCUCAGGUUCAAGUGGUUCUCCUUCCUCAGCCUCCC"
                "AAGUAGCUGGGAUUACAGGCACCUACCACCAGGGCCAGCUAAUUUUUGUA"
                "UGUUUAGUAGAAACGGGGUUUCACCAUGUUGGCCAGGCUGUUCUCGAACU"
                "CCUGACCUCAGGUGAUCUGCCCGCCUCGGUCUCCCAAAGUGCUGGGAUUA"
                "CAGGCAUGAGCCACUGCACCUGACCUGCUGAAUUGUUUAUAAUGGCAAGA"
                "AAUAGGAAACCCCCCAAUGUCUGUUGAACAGCUAUCACGUUGAACCACGU"
                "GAAACUGCUGUUUUCUAGGCCAAAAAUGGUGAGCGAUCAUUUAUUUCAUG"
                "AUUCAACCUGAUACAUUUACAUAGUGCAAAACUGUGUCACAGUUUCAGGC"
                "UUUUAUGAGGAAAGCGUUUCUGUGUAGAAACUGGAAGCUGUUCAGGGCAU"
                "CGGCAGCUGAACCCUGCUCCGUUGGUCAGCGUUACUAUCAUCUCGGAUCA"
                "UAUGGAGCUCAUGUCAGCCGUGUGGGUGGCGGGUGCACAGAGACGGUCUG"
                "GAAGGAAACACGCGGAUCUGAACAGCAGUAAUCCUGGGGGAUACGGGGGU"
                "UGGGCUAGAUUACAGAGGGCUCAUUUUCUACGUCAUGUAUUUUAUGAUAC"
                "UUGAAUUUUUUGAAAUGGGCAUUUAUUUUAUAACAUGUUAAAAUGUACUU"
                "UUUAAAUUAAGUCAUUUUGUAAUAUUUGAAUUUUUACAUUUGUUGUACAA"
                "UCAGGAAAAGCAAUAAAGAUUUUUCAAAAAUAG";
        EXPECT_STREQ(seq5, seqan::toCString((seqan::CharString)mrna_seqs[4]));

        const char *seq6 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCCCAAUGUCUUCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                "CCCCCCCUGGAAAAAACCUUAAAAAAUGUUUCCUCCAAAUCUGAUUUCAU"
                "UACAUUUCUGAAUUGUUGGGGUUUUUUUUGUUGUUUUGUUUUGUUUUGUA"
                "GAUGGAGUUUCACUUUUGUUGCCCAGGCUGGAGUGUAGUGGCGCGAUCUC"
                "GGCUCAGCCUCCCGAGUAGCUGGGAUUACAGGCAUGUGCCACCACGCCCG"
                "GCUAAUGUUUGUAUUUUUAGUAGAGACGGGGUUUCACCAUGUUGGUCAGG"
                "CUGGUCUCAAACUCCUGACCUCAGGUGAUCCGCCCACCUCAGCCUCCCAA"
                "AGUGCUGGGAUGACAGGUGUGAGCCACUGCGCCCAGCCUGAAUCAUUUCU"
                "UAUACCUUCUGACAGCCCAACUUCCAGAGGACAGCUCUGGGGUACUCGUU"
                "GGAUGUCUGUGAGUACCUGGUCAUACGGGUCAGUAGGGAUAAGAAUUGUC"
                "UCUGGGCUGAGGAAUUCUUCUGUUCUCUGGUUUCACCAGCGUUGGGUUUG"
                "CUCAUGUAAUGUGGUCACCAUACUCAAAUGGUGUCAUGGCUGAAGUUGGC"
                "CACCUUGCUUGAGGGACAAGUUGUUUAUGUAUCAGCUCUCUGCUGGGUCU"
                "CCCUUUCCAUGGCAAAUGGGCAGCUCCAUCCUCUUGACUCUUCUAAAUGC"
                "CCAAAAGAGGUGUCAUGCUUUGGGGGUACGAUGUUUAUACUCCGUAAAGA"
                "ACAUACAAGGACAUUCACUGCUGAUUUUUUUUUUUGUUUGUUUGAGACAG"
                "GGUCUCACUCUGUCGCUCAGGCUGGAGUGCAGUGAUGCAAUCUUGGCUCA"
                "CUGCAACCUCCGCCUCUCAGGUUCAAGUGGUUCUCCUUCCUCAGCCUCCC"
                "AAGUAGCUGGGAUUACAGGCACCUACCACCAGGGCCAGCUAAUUUUUGUA"
                "UGUUUAGUAGAAACGGGGUUUCACCAUGUUGGCCAGGCUGUUCUCGAACU"
                "CCUGACCUCAGGUGAUCUGCCCGCCUCGGUCUCCCAAAGUGCUGGGAUUA"
                "CAGGCAUGAGCCACUGCACCUGACCUGCUGAAUUGUUUAUAAUGGCAAGA"
                "AAUAGGAAACCCCCCAAUGUCUGUUGAACAGCUAUCACGUUGAACCACGU"
                "GAAACUGCUGUUUUCUAGGCCAAAAAUGGUGAGCGAUCAUUUAUUUCAUG"
                "AUUCAACCUGAUACAUUUACAUAGUGCAAAACUGUGUCACAGUUUCAGGC"
                "UUUUAUGAGGAAAGCGUUUCUGUGUAGAAACUGGAAGCUGUUCAGGGCAU"
                "CGGCAGCUGAACCCUGCUCCGUUGGUCAGCGUUACUAUCAUCUCGGAUCA"
                "UAUGGAGCUCAUGUCAGCCGUGUGGGUGGCGGGUGCACAGAGACGGUCUG"
                "GAAGGAAACACGCGGAUCUGAACAGCAGUAAUCCUGGGGGAUACGGGGGU"
                "UGGGCUAGAUUACAGAGGGCUCAUUUUCUACGUCAUGUAUUUUAUGAUAC"
                "UUGAAUUUUUUGAAAUGGGCAUUUAUUUUAUAACAUGUUAAAAUGUACUU"
                "UUUAAAUUAAGUCAUUUUGUAAUAUUUGAAUUUUUACAUUUGUUGUACAA"
                "UCAGGAAAAGCAAUAAAGAUUUUUCAAAAAUAG";
        EXPECT_STREQ(seq6, seqan::toCString((seqan::CharString)mrna_seqs[5]));
    }
}
