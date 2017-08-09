#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_site_score.hpp"   // PITASiteScores

using namespace seqan;

namespace ptddg {

//
// PITADGDuplexScores methods
//
void PITADGDuplexScores::clear_scores() {
    clear(mEffectiveSites);
}

int PITADGDuplexScores::calc_scores(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites) {

    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    String<int> const &mismatchPos = pSeedSites.get_mismatched_pos();
    std::string inputMiRNASeq;
    std::string inputMRNASeq;
    std::vector<int> inputMatchSeq;
    int seqStart = 0;
    int seqEnd = 0;

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    mVRws.preppare_duplexfold((int) length(pSeedSites.mEffectiveSites));

    inputMiRNASeq.resize(length(pMiRNASeq));
    create_input_mirna_seq(pMiRNASeq, inputMiRNASeq);

    inputMRNASeq.resize(TARGET_SEQ_LEN);

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            mEffectiveSites[i] = false;
            continue;
        }

        seqEnd = sitePos[i] + (mikan::SEEDLEN + 1);
        seqStart = seqEnd - TARGET_SEQ_LEN;
        if (seqStart < 0) {
            seqStart = 0;
        }

        create_input_mrna_seq(pMRNASeqs[mRNAPos[i]], seqStart, seqEnd, inputMRNASeq);
        create_input_matched_seq(seedTypes[i], mismatchPos[i], inputMatchSeq);
//        print_input(seedTypes[i], inputMiRNASeq, inputMRNASeq, inputMatchSeq);

        mVRws.duplexfold(i, inputMiRNASeq, inputMRNASeq, inputMatchSeq, inputMatchSeq);
//        mVRws.print_duplexfold_ret_vals(i);

        mAlign.create_align(i, pMiRNASeq, pMRNASeqs[mRNAPos[i]], seedTypes[i], (unsigned) sitePos[i],
                            mismatchPos[i]);

        mEffectiveSites[i] = true;
    }

    return 0;
}

void PITADGDuplexScores::create_input_mirna_seq(
        mikan::TRNAStr const &pMiRNASeq,
        std::string &pInputMiRNASeq) {
    for (unsigned i = 0; i < length(pMiRNASeq); ++i) {
        pInputMiRNASeq[i] = pMiRNASeq[i];
    }
}

void PITADGDuplexScores::create_input_mrna_seq(
        mikan::TRNAStr const &pMRNASeq,
        int pStart,
        int pEnd,
        std::string &pInputMRNASeq) {
    int seqStartPos = TARGET_SEQ_LEN - (pEnd - pStart);

    for (int i = 0; i < seqStartPos; ++i) {
        pInputMRNASeq[i] = 'A';
    }

    for (int i = 0; i < (pEnd - pStart); ++i) {
        pInputMRNASeq[seqStartPos + i] = pMRNASeq[pStart + i];
    }

}

void PITADGDuplexScores::create_input_matched_seq(
        CharString const &pSeedType,
        int pMismatchPos,
        std::vector<int> &pInputMatchSeq) {
    unsigned seedLen = lexicalCast<unsigned>(pSeedType[0]);
    pInputMatchSeq.resize((seedLen + 1));

    for (int i = 0; i < (int) seedLen + 1; ++i) {
        if (i == 0) {
            pInputMatchSeq[i] = 0;
        } else {
            pInputMatchSeq[i] = 1 + i;
        }
    }

    if (pSeedType == "8mer_MM" || pSeedType == "7mer_MM" || pSeedType == "8mer_MMGU" || pSeedType == "7mer_MMGU") {
        pInputMatchSeq[mikan::SEEDLEN - pMismatchPos] = 0;
    }

}

void PITADGDuplexScores::print_input(
        CharString const &pSeedType,
        std::string &pInputMiRNASeq,
        std::string &pInputMRNASeq,
        std::vector<int> &pInputMatchSeq) {
    std::cout << "seed type:   " << toCString(pSeedType) << std::endl;
    std::cout << "miRNA seq:   " << pInputMiRNASeq << std::endl;
    std::cout << "mRNA seq:    " << pInputMRNASeq << std::endl;
    std::cout << "constraints: ";
    for (unsigned k = 0; k < pInputMatchSeq.size(); ++k) {
        std::cout << pInputMatchSeq[k];
        if (k != pInputMatchSeq.size() - 1) {
            std::cout << ",";
        } else {
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

//
// PITADGOpenScores methods
//

void PITADGOpenScores::clear_scores() {
    clear(mEffectiveSites);
}

int PITADGOpenScores::calc_scores(
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites) {

    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    std::string inputMRNASeq;
    int seqStart = 0;
    int seqEnd = 0;
    int seqStart2, seqEnd2;
    int paramU, paramS, paramFT;
    int seqLen;

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));

    paramU = mFlankUp;
    paramS = DDG_OPEN + DDG_AREA + mFlankDown - 1;
    paramFT = DDG_OPEN + mFlankDown;
    mVRws.prepare_ddg4((int) length(pSeedSites.mEffectiveSites), paramU, paramS, paramFT, paramFT);

    seqLen = DDG_OPEN + DDG_AREA * 2 + mFlankUp + mFlankDown;
    resize(inputMRNASeq, seqLen);

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            mEffectiveSites[i] = false;
            continue;
        }

        seqEnd = sitePos[i] + (mikan::SEEDLEN + 1);
        seqStart = seqEnd - DDG_OPEN;
        seqStart2 = seqStart - (DDG_AREA + mFlankDown);
        seqEnd2 = seqEnd + (DDG_AREA + mFlankUp);

        create_input_mrna_seq(pMRNASeqs[mRNAPos[i]], seqStart2, seqEnd2, inputMRNASeq);
//        print_input(inputMRNASeq);

        mVRws.calc_ddg4(i, inputMRNASeq);
//        mVRws.print_ddg4_ret_vals(i);

        mEffectiveSites[i] = true;

    }

    return 0;
}

void PITADGOpenScores::create_input_mrna_seq(
        mikan::TRNAStr const &pMRNASeq,
        int pStart,
        int pEnd,
        std::string &pInputMRNASeq) {
    int k = 0;
    for (int i = pStart; i < pEnd; ++i) {
        if ((i < 0) || (i >= (int) length(pMRNASeq))) {
            pInputMRNASeq[k] = 'A';
        } else {
            pInputMRNASeq[k] = pMRNASeq[i];
        }
        ++k;
    }
}

void PITADGOpenScores::print_input(std::string &pInputMRNASeq) {
    std::cout << "mRNA seq: " << pInputMRNASeq << std::endl;
}

//
// PITASiteScores methods
//
void PITASiteScores::init_from_args() {
    set_backtrack(mOpts.mOutputAlign);

    mDGOpenScores.mFlankUp = mOpts.mFlankUp;
    mDGOpenScores.mFlankDown = mOpts.mFlankDown;
}

void PITASiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    clear(mDDGScores);
    mDGDuplexScores.clear_scores();
    mDGOpenScores.clear_scores();
    mAlign.clear_align();
}

int PITASiteScores::calc_scores(
        mikan::TRNAStr const &miRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &) {

    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    resize(mDDGScores, length(pSeedSites.mEffectiveSites));
    mAlign.resize_align((int) length(pSeedSites.mEffectiveSites));

    mDGDuplexScores.calc_scores(miRNASeq, pMRNASeqs, pSeedSites);
    mDGOpenScores.calc_scores(pMRNASeqs, pSeedSites);

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i] || !mDGDuplexScores.mEffectiveSites[i] ||
            !mDGOpenScores.mEffectiveSites[i]) {
            mEffectiveSites[i] = false;
            pSeedSites.mEffectiveSites[i] = false;
            mDDGScores[i] = 0.0;
        } else {
            mDDGScores[i] = mVRws.get_dgall(i) + (mVRws.get_dg1(i) - mVRws.get_dg0(i));
            mEffectiveSites[i] = true;
        }

    }

    return 0;
}

void PITASiteScores::print_alignment(int pIdx) {
    std::stringstream stream;

    stream << "mRNA   5' " << mAlign.mAlignMRNA[pIdx] << " 3'";
    stream << std::endl;
    stream << "          " << mAlign.mAlignBars[pIdx] << "   ";
    stream << std::endl;
    stream << "miRNA  3' " << mAlign.mAlignMiRNA[pIdx] << " 5'";
    stream << std::endl;


    std::cout << stream.str();
}

} // namespace ptddg
