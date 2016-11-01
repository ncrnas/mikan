#include <iostream>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class U3009 : public TestIO
    {
    protected:
        U3009() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_009.fasta";
            OFNAME1 = (char *)"test_output_1.txt";
            OFNAME2 = (char *)"test_ds1.txt";
            OMPATH = (char *)"mkmiranda/";
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

    TEST_F(U3009, mrna_fasta) {
        read_files();
        mrna_ids = coreInput.get_mrna_ids();
        mrna_seqs = coreInput.get_mrna_seqs();

        EXPECT_EQ(1u, length(mrna_ids));

        const char *id1 = "hg18_refGene NM_001259 range=chr7:92072173-92082389 5'pad=0 3'pad=0 "
                "revComp=TRUE strand=- repeatMasking=none";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *seq1 = "GGCCUCAGCAGCCGCCUUAAGCUGAUCCUGCGGAGAACACCCUUGGUGGC"
                "UUAUGGGUCCCCCUCAGCAAGCCCUACAGAGCUGUGGAGGAUUGCUAUCU"
                "GGAGGCCUUCCAGCUGCUGUCUUCUGGACAGGCUCUGCUUCUCCAAGGAA"
                "ACCGCCUAGUUUACUGUUUUGAAAUCAAUGCAAGAGUGAUUGCAGCUUUA"
                "UGUUCAUUUGUUUGUUUGUUUGUCUGUUUGUUUCAAGAACCUGGAAAAAU"
                "UCCAGAAGAAGAGAAGCUGCUGACCAAUUGUGCUGCCAUUUGAUUUUUCU"
                "AACCUUGAAUGCUGCCAGUGUGGAGUGGGUAAUCCAGGCACAGCUGAGUU"
                "AUGAUGUAAUCUCUCUGCAGCUGCCGGGCCUGAUUUGGUACUUUUGAGUG"
                "UGUGUGUGCAUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUGUAUGUGAG"
                "AGAUUCUGUGAUCUUUUAAAGUGUUACUUUUUGUAAACGACAAGAAUAAU"
                "UCAAUUUUAAAGACUCAAGGUGGUCAGUAAAUAACAGGCAUUUGUUCACU"
                "GAAGGUGAUUCACCAAAAUAGUCUUCUCAAAUUAGAAAGUUAACCCCAUG"
                "UCCUCAGCAUUUCUUUUCUGGCCAAAAGCAGUAAAUUUGCUAGCAGUAAA"
                "AGAUGAAGUUUUAUACACACAGCAAAAAGGAGAAAAAAUUCUAGUAUAUU"
                "UUAAGAGAUGUGCAUGCAUUCUAUUUAGUCUUCAGAAUGCUGAAUUUACU"
                "UGUUGUAAGUCUAUUUUAACCUUCUGUAUGACAUCAUGCUUUAUCAUUUC"
                "UUUUGGAAAAUAGCCUGUAAGCUUUUUAUUACUUGCUAUAGGUUUAGGGA"
                "GUGUACCUCAGAUAGAUUUUAAAAAAAAGAAUAGAAAGCCUUUAUUUCCU"
                "GGUUUGAAAUUCCUUUCUUCCCUUUUUUUGUUGUUGUUAUUGUUGUUUGU"
                "UGUUGUUAUUUUGUUUUUGUUUUUAGGAAUUUGUCAGAAACUCUUUCCUG"
                "UUUUGGUUUGGAGAGUAGUUCUCUCUAACUAGAGACAGGAGUGGCCUUGA"
                "AAUUUUCCUCAUCUAUUACACUGUACUUUCUGCCACACACUGCCUUGUUG"
                "GCAAAGUAUCCAUCUUGUCUAUCUCCCGGCACUUCUGAAAUAUAUUGCUA"
                "CCAUUGUAUAACUAAUAACAGAUUGCUUAAGCUGUUCCCAUGCACCACCU"
                "GUUUGCUUGCUUUCAAUGAACCUUUCAUAAAUUCGCAGUCUCAGCUUAUG"
                "GUUUAUGGCCUCGAUUCUGCAAACCUAACAGGGUCACAUAUGUUCUCUAA"
                "UGCAGUCCUUCUACCUGGUGUUUACUUUUGUUACCUAAAUAAUGAGUAGG"
                "AUCUUGUUUUGUUUUAUCACCAGCACACAGAUUGCUAUAAACUGUUACUU"
                "UGUGAAUUACAUUUUUAUAGAAGAUAUUUUCAGUGUCUUUACCUGAGGGU"
                "AUGUCUUUAGCUAUGUUUUAGGGCCAUACAUUUACUCUAUCAAAUGAUCU"
                "UUUCUCCAUCCCCCAGGCUGUGCUUAUUUCUAGUGCCUUGUGCUCACUCC"
                "UGCUCUCUACAGAGCCAGCCUGGCCUGGGCAUUGUAAACAGCUUUUCCUU"
                "UUUCUCUUACUGUUUUCUCUACAGUCCUUUAUAUUUCAUACCAUCUCUGC"
                "CUUAUAAGUGGUUUAGUGCUCAGUUGGCUCUAGUAACCAGAGGACACAGA"
                "AAGUAUCUUUUGGAAAGUUUAGCCACCUGUGCUUUCUGACUCAGAGUGCA"
                "UGCAACAGUUAGAUCAUGCAACAGUUAGAUUAUGUUUAGGGUUAGGAUUU"
                "UCAAAGAAUGGAGGUUGCUGCACUCAGAAAAUAAUUCAGAUCAUGUUUAU"
                "GCAUUAUUAAGUUGUACUGAAUUCUUUGCAGCUUAAUGUGAUAUAUGACU"
                "AUCUUGAACAAGAGAAAAAACUAGGAGAUGUUUCUCCUGAAGAGCUUUUG"
                "GGGUUGGGAACUAUUCUUUUUUAAUUGCUGUACUACUUAACAUUGUUCUA"
                "AUUCAGUAGCUUGAGGAACAGGAACAUUGUUUUCUAGAGCAAGAUAAUAA"
                "AGGAGAUGGGCCAUACAAAUGUUUUCUACUUUCGUUGUGACAACAUUGAU"
                "UAGGUGUUGUCAGUACUAUAAAUGCUUGAGAUAUAAUGAAUCCACAGCAU"
                "UCAAGGUCAGGUCUACUCAAAGUCUCACAUGGAAAAGUGAGUUCUGCCUU"
                "UCCUUUGAUCGAGGGUCAAAAUACAAAGACAUUUUUGCUAGGGCCUACAA"
                "AUUGAAUUUAAAAACUCACUGCACUGAUUCAUCUGAGCUUUUUGGUUAGU"
                "AUUCAUGGCUAGAGUGAACAUAGCUUUAGUUUUUGCUGUUGUAAAAGUGU"
                "UUUCAUAAGUUCACUCAAGAAAAAUGCAGCUGUUCUGAACUGGAAUUUUU"
                "CAGCAUUCUUUAGAAUUUUAAAUGAGUAGAGAGCUCAACUUUUAUUCCUA"
                "GCAUCUGCUUUUGACUCAUUUCUAGGCAGUGCUUAUGAAGAAAAAUUAAA"
                "GCACAAACAUUCUGGCAUUCAAUCGUUGGCAGAUUAUCUUCUGAUGACAC"
                "AGAAUGAAAGGGCAUCUCAGCCUCUCUGAACUUUGUAAAAAUCUGUCCCC"
                "AGUUCUUCCAUCGGUGUAGUUGUUGCAUUUGAGUGAAUACUCUCUUGAUU"
                "UAUGUAUUUUAUGUCCAGAUUCGCCAUUUCUGAAAUCCAGAUCCAACACA"
                "AGCAGUCUUGCCGUUAGGGCAUUUUGAAGCAGAUAGUAGAGUAAGAACUU"
                "AGUGACUACAGCUUAUUCUUCUGUAACAUAUGGUUUCAAACAUCUUUGCC"
                "AAAAGCUAAGCAGUGGUGAACUGAAAAGGGCAUAUUGCCCCAAGGUUACA"
                "CUGAAGCAGCUCAUAGCAAGUUAAAAUAUUGUGACAGAUUUGAAAUCAUG"
                "UUUGAAUUUCAUAGUAGGACCAGUACAAGAAUGUCCCUGCUAGUUUCUGU"
                "UUGAUGUUUGGUUCUGGCGGCUCAGGCAUUUUGGGAACUGUUGCACAGGG"
                "UGGAGUCAAAACAACCUACAUAUAAAAAAGAAACUUGUCCAUUUAGCUUU"
                "CAUAAGAAAUCCCAUGGCAAAGGGUAAUAAAAAGGACCUAAUCUUAAAAA"
                "UACAAUUUCUAAGCACUUGUAAGAACCCAGUGGGUUGGAGCCUCCCACUU"
                "UGUCCCUCCUUUGAAGUGGAUGGGAACUCAAGGUGCAAAGAACCUGUUUU"
                "GGAAGAAAGCUUGGGGCCAUUUCAGCCCCCUGUAUUCUCAUGAUUUUCUC"
                "UCAGGAAGCACACACUGUGAAUGGCAGACUUUUCAUUUAGCCCCAGGUGA"
                "CUUACUAAAAAUAGUUGAAAAUUAUUCACCUAAGAAUAGAAUCUCAGCAU"
                "UGUGUUAAAUAAAAAUGAAAGCUUUAGAAGGCAUGAGAUGUUCCUAUCUU"
                "AAAUAAAGCAUGUUUCUUUUCUAUAGAGAAAUGUAUAGUUUGACUCUCCA"
                "GAAUGUACUAUCCAUCUUGAUGAGAAAACUCUUAAAUAGUACCAAACAUU"
                "UUGAACUUUAAAUUAUGUAUUUAAAGUGAGUGUUUAAGAAACUGUAGCUG"
                "CUUCUUUUACAAGUGGUGCCUAUUAAAGUCAGUAAUGGCCAUUAUUGUUC"
                "CAUUGUGGAAAUUAAAUUAUGUAAGCUUCCUAAUAUCAUAAACAUAUUAA"
                "AAUUCUUCUAAAAUAUUGCUUUUCUUUUAAGUGACAAUUUGACUAUUCUU"
                "AUGAUAAGCACAUGAGAGUGUCUUACAUUUUCCAAAAGCAGGCUUUAAUU"
                "GCAUAGUUGAGUCUAGGAAAAAAUAAUGUUAAAAGUGAAUAUGCCACCAU"
                "AAUUACUUAAUUAUGUUAGUAUAGAAACUACAGAAUAUUUACCCUGGAAA"
                "GAAAAUAUUGGAAUGUUAUUAUAAACUCUUAGAUAUUUAUAUAAUUCAAA"
                "AGAAUGCAUGUUUCACAUUGUGACAGAUAAAGAUGUAUGAUUUCUAAGGC"
                "UUUAAAAAUUAUUCAUAAAACAGUGGGCAAUAGAUAAAGGAAAUUCUGGA"
                "GAAAAUGAAGGUAUUUAAAGGGUAGUUUCAAAGCUAUAUAUAUUUUGAAG"
                "GAUAUAUUCUUUAUGAACAAAUAUAUUGUAAAAAUUUAUACUAAGGUCAU"
                "CUGGUAACUGUGGGAUUAAUAUGGUCGAAAACAAAUGUUAUGGAGAAGCU"
                "GUCCCAAGCAAACUAAAUUACCUGUACUUUUUUCCCAUUUCAAGGGAAGA"
                "GGCAACCACAUGAAGCAAUACUUCUUACACAUGCCUAAGAACGUUCAUUG"
                "AAAAAAUAAAUUUUUAAAAGGCAUGUGUUUCCUAUGCCACCAAUACUUUU"
                "GAAAAAUUGUGAACCUUACCCAAAACCAUUUAUCAUGUCCAUUAAGUAUA"
                "UUUGGGUAUAUAAUUAGGAAGAUAUUUACAUGUUCCAUCUCCACAGUGGA"
                "AAAACUUAUUGAGGCUACCAAAGUGUGCCAAGAAAUGUAAGUCCUUAGAG"
                "UAAUUAGAAAUGCUGUUUUCCUCAAAAGCAUGAGAAACUAGCAUUUUCAU"
                "UUCUUAUUUACUCCCUUUCUAUAUCAAUGCAAUUCACAACCCAAUUUUAA"
                "UACAUCCCUAUAUCUCAAGCAUUUCUAUCUUGUACUUUUUCAGAAAAUAA"
                "ACCAAAAAUAAUCCUUUGGUCUCUCUAUCUUCUGACCUUUGUAAGCAACA"
                "GAAAUGUAAAAACAGAAGGGGUCCAAUUUUUACACGUUUUUUUCUCAAGU"
                "AGCCUUUCUGGGGAUUUUUAUUUUCUUAAUGAAGUGCCAAUCAGCUUUUC"
                "AAAAUGUUUUCUAUUUCUCAGCAUUUCCAGGAAGUGAUAACGUUUAGCUA"
                "AAUGAGUAGAAGUGGACUUCCUUCAACAUAUUGUUACCUUGUCUAGCCUU"
                "AGGAAGAAAACAAGAGCCACCUGAAAAUAAAUACAGGCUCUUUUCGAGCA"
                "UCUGCUGAAAUACUGUUACAGCAAUUUGAAGUUGAUGUGGUAGGAAAGGA"
                "AGGUGACUUUUCUUGCAAAAGUCUUUCUAAACAUUCACACUGUCCUAAGA"
                "GAUGAGCUUUCUUGUUUUAUUCCGGUAUAUUCCACAAGGUGGCACUUUUA"
                "GAGAAAAACAAAUCUGAUGAAGACUAAAGAGGUACUUCUAAAAGAGAUUU"
                "CAUUCUAACUUUAUUUUUCUGCGCAUAUUUAACUCUUUCCUAGCACUUGU"
                "UUUUUGGGAUGAUUAAUAGUCUCUAUAAUGUUCUGUAACUUCAAUAUUUU"
                "ACUUGUUACCUAGGUUCUGAACAAUUGUCUGCAAAUAAAUUGUUCUUAAG"
                "GAUGGAUAAUACACCCAUUUUGAUCAUUUAAGUAAAGAAAGCCUAGUCAU"
                "UCAUUCAGUCAAGAAAAAAUUUUUGAAGUACCCAGUUACCUUACUUUUCU"
                "AGAUUAAAACAGGCUUAGUUACUAAAAAGGCAGUCCUCAUCUGUGAACAG"
                "GAUAGUUUCGUUAGAAGUAUAAAACUCCUUUAGUGGCCCCAGUUAAAACA"
                "CACAUACCCUCUCUGCUGCUUUCAAAUUCCCUAGCAUGGUGGCCUUUCAA"
                "CAUUGAUUAAAUUUUAAAAUCCUAAUUUAAAGAUCAGGUGAGCAAAAUGA"
                "GUAGCACAUCAGUAAUUCAGUAGACAAAACUUUUGUCUGAAAAAUUGCUG"
                "UAUUGAAACAGAGCCCUAAAAUACCAAAAGACCAGGUAAUUUUAACAUUU"
                "GUGGAAUCACAAAUGUAAAUUCAUAAGAAGCUCUAAUUAAAAAAAAAAAG"
                "UCUGAAGUAUAUGAGCAUAACAACUUAGGAGUGUGUCUACAUACUUAACU"
                "UUUGAAGUUUUUUGGCAACUUUAUAUACUUUUUUUAAAUUUACAAGUCUA"
                "CUUAAAGACUUCUUAUACCCCAAAUGAUUAAGUUAAUUUUAGAGGUCACC"
                "UUUCUCACAGCAGUGUCACUUGAAAUUUAGUAGGGAAGGAUAUUGCAGUA"
                "UUUUUCAGUUUCCUUAGCACAGCACCACAGAAAGCAGCUUAUUCCUUUUG"
                "AGUGGCAGACACUCGACGGUGCCUGCCCAACUUUCCUCCUGAGUGGCAAG"
                "CAGAUGAGUCUCAGUAAUUCAUACUGAACCAAAAUGCCACAUACACUAGG"
                "GGCAGUCAGAAACUGGCUGAGAAAUCCCCCGCCUCAUUCGCCCCUCUGCU"
                "CCCAGGAACUAGAGUCCAGUUAAAGCCCCUAUGCGAAAGGCCGAAUUCCA"
                "CCCCAGGGUUUGUUAUAACAGUGGCCAGUCUGAACCCCAUUUGCUCGUGC"
                "UCAAAACUUGAUUCCCACUUGAAAGCCUUCCGGGCGCGCUGCCUCGUUGC"
                "CCCGCCCCUUUGGCAGGAGAGAGGCAGUGGGCGAGGCCGGGCUGGGGCCC"
                "CGCCUCCCACUCACCUGCCGGUGCCUGAAAUUAUGUGCGGCCCCGCGGGC"
                "UGCUUUCCGAGGUCAGAGUGCCCUGCUGCUGUCUCAGAGGCAUCUGUUCU";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }
}