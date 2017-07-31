#include "mk_typedef.hpp"   // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_score.hpp"   // PITASiteScores, PITATotalScores

using namespace seqan;

namespace ptddg {

//
// PITAAlign methods
//
void PITAAlign::clear_align() {
    clear(mEffectiveSites);
    clear(mAlignMRNA);
    clear(mAlignBars);
    clear(mAlignMiRNA);
}

void PITAAlign::resize_align(unsigned pSize) {
    resize(mEffectiveSites, pSize, false);
    resize(mAlignMRNA, pSize);
    resize(mAlignBars, pSize);
    resize(mAlignMiRNA, pSize);
}

void PITAAlign::create_align(
        int pId,
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNAStr const &pMRNASeq,
        CharString const &pSeedType,
        unsigned pSitePos,
        int) {
    int seedLen = lexicalCast<int>(pSeedType[0]);
    int startPos;
    int pos;
    char mChar;


    resize(mAlignMiRNA[pId], length(pMiRNASeq));
    for (unsigned i = 0; i < length(pMiRNASeq); ++i) {
        mAlignMiRNA[pId][length(pMiRNASeq) - i - 1] = pMiRNASeq[i];
    }

    resize(mAlignMRNA[pId], length(pMiRNASeq), ' ');
    startPos = pSitePos + (INDEXED_SEQ_LEN + 1) - (int) length(pMiRNASeq);
    for (unsigned i = 0; i < length(pMiRNASeq); ++i) {
        if (startPos + (int) i > 0) {
            mAlignMRNA[pId][i] = pMRNASeq[startPos + i];
        }
    }

    resize(mAlignBars[pId], length(pMiRNASeq), ' ');
    for (int i = 1; i < seedLen + 1; ++i) {
        pos = (int) length(pMiRNASeq) - 1 - i;
        mChar = ' ';
        if ((mAlignMiRNA[pId][pos] == 'C' && mAlignMRNA[pId][pos] == 'G')
            || (mAlignMiRNA[pId][pos] == 'G' && mAlignMRNA[pId][pos] == 'C')
            || (mAlignMiRNA[pId][pos] == 'A' && mAlignMRNA[pId][pos] == 'U')
            || (mAlignMiRNA[pId][pos] == 'U' && mAlignMRNA[pId][pos] == 'A')) {
            mChar = '|';
        } else if ((mAlignMiRNA[pId][pos] == 'G' && mAlignMRNA[pId][pos] == 'U')
                   || (mAlignMiRNA[pId][pos] == 'U' && mAlignMRNA[pId][pos] == 'G')) {
            mChar = ':';
        }

        mAlignBars[pId][pos] = mChar;
    }

}

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

        seqEnd = sitePos[i] + (INDEXED_SEQ_LEN + 1);
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
        pInputMatchSeq[INDEXED_SEQ_LEN - pMismatchPos] = 0;
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

        seqEnd = sitePos[i] + (INDEXED_SEQ_LEN + 1);
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
        mikan::MKSeedSites &pSeedSites) {

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

//
// PITATotalScores methods
//
void PITATotalScores::clear_scores() {
    clear(mMRNAPos);
    clear(mSiteNum);
    clear(mTotalScores);
}

int PITATotalScores::calc_scores(
        PITASeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        PITASiteScores &pSiteScores) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mTotalScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs));

    float score, max_score, exp_diff, total_score;
    unsigned site_count, max_idx;
    mikan::TSitePosSet sitePos;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        score = 0;
        max_score = -FLT_MAX;
        site_count = 0;
        max_idx = 0;

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }
            appendValue(sitePos, rnaSitePosMap[i][j]);
            score = -1.0 * pSiteScores.get_score(rnaSitePosMap[i][j]);
            if (score > max_score) {
                max_score = score;
                max_idx = j;
            }
            ++site_count;
        }

        if (length(sitePos) == 0) {
            continue;
        }

        total_score = 1.0;
        for (unsigned j = 0; j < length(sitePos); ++j) {
            if (j == max_idx) {
                continue;
            }
            score = -1.0 * pSiteScores.get_score(sitePos[j]);
            exp_diff = score - max_score;
            if (exp_diff > MIN_EXP_DIFF) {
                total_score += std::exp(exp_diff);
            }
        }

        mTotalScores[i] = -1.0 * (max_score + std::log(total_score));
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = site_count;

        clear(sitePos);
    }

    return 0;

}

} // namespace ptddg
