#include "mk_typedef.hpp"  // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_score.hpp"   // RH2SiteScores, RH2TotalScores

using namespace seqan;

namespace rh2mfe {

//
// RH2SiteScores methods
//
void RH2SiteScores::init_from_args() {
    mOverlapDef = mOpts.mOverlapDef;
}

void RH2SiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    clear(mMFEScores);
    clear(mNormScores);
    mRHRetVals.clear();
}

int RH2SiteScores::calc_scores(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites) {
    
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    mikan::TRNAStr miRNASeq = pMiRNASeq;
    mikan::TRNAStr mRNASeq;
    std::vector<char> rhMiRNASeq;
    std::vector<char> rhMRNASeq;
    int seqStart = 0;
    int seqEnd = 0;
    int targetLen, queryLen;

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    resize(mMFEScores, length(pSeedSites.mEffectiveSites));
    resize(mNormScores, length(pSeedSites.mEffectiveSites));
    mRHRetVals.resize(length(pSeedSites.mEffectiveSites));

    reverse(miRNASeq);
    queryLen = length(miRNASeq);
    rhMiRNASeq.resize((unsigned) (queryLen + 2));
    create_rh_seq(miRNASeq, rhMiRNASeq);
    rhMRNASeq.resize((unsigned) (mMaxMRNALen + 2));

    mRHCore.init_workspace(mMaxMRNALen, queryLen);

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            mEffectiveSites[i] = false;
            mMFEScores[i] = 0.0;
            mNormScores[i] = 0.0;
            continue;
        }

        seqStart = sitePos[i] - (mMaxMRNALen - 8);
        seqEnd = sitePos[i] + 8;

        if (length(pMRNASeqs[mRNAPos[i]]) < (unsigned) seqEnd) {
            seqEnd = (int) length(pMRNASeqs[mRNAPos[i]]);
        }

        if (seqStart < 0) {
            seqStart = 0;
        }

        mRNASeq = infix(pMRNASeqs[mRNAPos[i]], seqStart, seqEnd);
        create_rh_seq(mRNASeq, rhMRNASeq);
        targetLen = length(mRNASeq);

//        write_seq_info(miRNASeq, mRNASeq, rhMiRNASeq, rhMRNASeq);

        mRHCore.set_target_seq(rhMRNASeq, targetLen);
        mRHCore.set_query_seq(rhMiRNASeq, queryLen);
        mRHCore.mainloop(mRHRetVals[i]);
        mRHRetVals[i].mTargetPos0 += seqStart;

//        write_alignment(i);

        mMFEScores[i] = mRHRetVals[i].mMfe;
        mEffectiveSites[i] = mRHRetVals[i].mEffective;
        if (!mEffectiveSites[i]) {
            pSeedSites.mEffectiveSites[i] = false;
        }
        calc_normalized_score(i, targetLen, queryLen);

    }

    return 0;
}

void RH2SiteScores::create_rh_seq(mikan::TRNAStr const &pRNASeq, std::vector<char> &pRHSeq) {
    pRHSeq[0] = 5;

    for (unsigned i = 0; i < length(pRNASeq); ++i) {
        pRHSeq[i + 1] = (char) ordValue(pRNASeq[i]);
    }

    pRHSeq[pRHSeq.size() - 1] = 0;

}

void RH2SiteScores::calc_normalized_score(int pIdx, int pTargetLen, int pQueryLen) {
    float mfx = mMFEScores[pIdx];
    mNormScores[pIdx] = -1.0f * mfx / log((float) (pTargetLen * pQueryLen));
}

void RH2SiteScores::write_alignment(int pIdx, bool align_only) {
    std::stringstream stream;

    if (!align_only) {
        stream << "position(target 5'): " << mRHRetVals[pIdx].mTargetPos0 + 1 << std::endl;
        stream << "# mfe: " << mRHRetVals[pIdx].mMfe << " kcal/mol" << std::endl;
        stream << "# normalized score: " << mNormScores[pIdx] << std::endl;
        stream << std::endl;
    }
    stream << "target 5' " << mRHRetVals[pIdx].mTargetSub.c_str() << " 3'" << std::endl;
    stream << "          " << mRHRetVals[pIdx].mTarget.c_str() << "   " << std::endl;
    stream << "          " << mRHRetVals[pIdx].mQeuery.c_str() << "   " << std::endl;
    stream << "miRNA  3' " << mRHRetVals[pIdx].mQeuarySub.c_str() << " 5'" << std::endl;

    std::cout << stream.str();
}

void RH2SiteScores::write_seq_info(
        mikan::TRNAStr &pMiSeq,
        mikan::TRNAStr &pMRNASeq,
        std::vector<char> &pRhMiRNASeq,
        std::vector<char> &pRhMRNASeq) {
    std::stringstream stream;

    stream << "miRNA: " << pMiSeq << "," << length(pMiSeq);
    stream << std::endl;
    stream << "miRNA: ";
    for (std::vector<char>::const_iterator j = pRhMiRNASeq.begin(); j != pRhMiRNASeq.end(); ++j) {
        stream << (int) *j;
    }
    stream << std::endl;

    stream << "mRNA:  " << pMRNASeq << "," << length(pMRNASeq);
    stream << std::endl;
    stream << "mRNA:  ";
    for (std::vector<char>::const_iterator j = pRhMRNASeq.begin(); j != pRhMRNASeq.end(); ++j) {
        stream << (int) *j;
    }
    stream << std::endl << std::endl;

    std::cout << stream.str();
}

//
// RH2TotalScores methods
//
void RH2TotalScores::clear_scores() {
    clear(mMRNAPos);
    clear(mSiteNum);
    clear(mTotalScores);
    clear(mTotalNormScores);
}

int RH2TotalScores::calc_scores(
        RH2SeedSites &pSeedSites,
        RH2SiteScores &pMFEScores,
        const seqan::String<unsigned> &pSortedSites) {

    const String<unsigned> &siteMRNAPos = pSeedSites.get_mrna_pos();
    int prevPos = -1;
    int newIdx = -1;
    String<int> newIdices;
    unsigned posIdx;

    resize(newIdices, length(siteMRNAPos));
    for (unsigned i = 0; i < length(pSortedSites); ++i) {
        posIdx = (unsigned) pSortedSites[i];

        if (!pMFEScores.mEffectiveSites[posIdx]) {
            continue;
        }

        if (prevPos != (int) siteMRNAPos[posIdx]) {
            ++newIdx;
        }
        newIdices[i] = newIdx;
        prevPos = (int) siteMRNAPos[posIdx];
    }

    resize(mTotalScores, newIdx + 1, 0.0);
    resize(mTotalNormScores, newIdx + 1, 0.0);
    resize(mMRNAPos, newIdx + 1);
    resize(mSiteNum, newIdx + 1, 0);

    for (unsigned i = 0; i < length(pSortedSites); ++i) {
        posIdx = (unsigned) pSortedSites[i];

        if (!pMFEScores.mEffectiveSites[posIdx]) {
            continue;
        }

        mMRNAPos[newIdices[i]] = siteMRNAPos[posIdx];
        mTotalScores[newIdices[i]] += pMFEScores.get_score(posIdx);
        mTotalNormScores[newIdices[i]] += pMFEScores.get_norm_score(posIdx);
        mSiteNum[newIdices[i]] += 1;
    }

    return 0;
}

} // namespace rh2mfe
