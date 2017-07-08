#include "mk_typedef.hpp"   // TRNATYPE
#include "mr3_score.hpp"    // MR3DDGScores, MR3TotalScores

using namespace seqan;

namespace mr3as {

//
// MR3AlignScores methods
//
void MR3AlignScores::clear_scores() {
    clear(mEffectiveSites);
    clear(mAlignScores);
}

int MR3AlignScores::calc_scores(
        MR3SeedSites &pSeedSites,
        mikan::TRNATYPE const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    const String<int> &mmPos = pSeedSites.get_mismatched_pos();
    const StringSet<CharString> &seedTypes = pSeedSites.get_seed_types();

    mikan::TRNATYPE iMiRNASeq;
    mikan::TRNATYPE iMRNASeq;
    mikan::TRNATYPE iMiRNASeedSeq;
    mikan::TRNATYPE iMRNASeedSeq;
    Rna5String iMiRNA3pSeq;
    Rna5String iMRNA3pSeq;
    int seqStart = 0;
    int seqEnd = 0;
    int score;
    int mm;
    bool noA1;

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    resize(mAlignScores, length(pSeedSites.mEffectiveSites));

    create_input_mirna_seq(pMiRNASeq, iMiRNASeq, iMiRNASeedSeq, iMiRNA3pSeq);

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
//        std::cout << seedTypes[i] << "," <<  sitePos[i] << ","<< pSeedSites.mEffectiveSites[i] << std::endl;
        if (!pSeedSites.mEffectiveSites[i]) {
            mEffectiveSites[i] = false;
            continue;
        }

        if (seedTypes[i] == "7mer_BT" || seedTypes[i] == "8mer_BT") {
            seqEnd = sitePos[i] + (INDEXED_SEQ_LEN + 2);
            seqStart = seqEnd - (TARGET_SEQ_LEN + 1);
            mm = mmPos[i];
        } else {
            seqEnd = sitePos[i] + (INDEXED_SEQ_LEN + 1);
            seqStart = seqEnd - TARGET_SEQ_LEN;
            mm = -1;
        }
        if (seqStart < 0) {
            seqStart = 0;
        }

        noA1 = false;
        create_input_mrna_seq(pMiRNASeq, pMRNASeqs[mRNAPos[i]], seqStart, seqEnd, seedTypes[i],
                              iMRNASeq, iMRNASeedSeq, iMRNA3pSeq, noA1);
//        print_input(iMiRNASeq, iMRNASeq);
//        print_input(iMiRNASeedSeq, iMRNASeedSeq);
//        std::cout << "miRNA seq:   " << length(iMiRNA3pSeq) << "," << iMiRNA3pSeq;
//        std::cout << std::endl;
//        std::cout << "mRNA seq2:    " << length(iMRNA3pSeq) << "," << iMRNA3pSeq;
//        std::cout << std::endl;

        mAlign.align_seed(i, iMiRNASeedSeq, iMRNASeedSeq, mm);

        mAlign.init_3p_align(i);
        if (length(iMRNA3pSeq) > 1) {
            mAlign.align_3p(i, iMiRNA3pSeq, iMRNA3pSeq);
        }

        score = mAlign.get_align_score(i);
//        std::cout << seedTypes[i] << "," <<  sitePos[i] << ","<< score << std::endl;

        if (mMinAlignScore > score) {
            mEffectiveSites[i] = false;
            pSeedSites.mEffectiveSites[i] = false;
        } else {
            mAlignScores[i] = (float) score;
            mAlign.combine_alignments(i, pMiRNASeq, iMRNASeq, noA1);
            mEffectiveSites[i] = true;
        }
    }

    return 0;
}

void MR3AlignScores::create_input_mirna_seq(
        mikan::TRNATYPE const &pMiRNASeq,
        mikan::TRNATYPE &pIMiRNASeq,
        mikan::TRNATYPE &pIMiRNASeedSeq,
        Rna5String &pIMiRNA3pSeq) {
    int idxseed = 0;
    int idx3p = 1;

    resize(pIMiRNASeq, length(pMiRNASeq));
    resize(pIMiRNASeedSeq, SEED_REGION_LEN - 1);
    resize(pIMiRNA3pSeq, length(pMiRNASeq) - SEED_REGION_LEN - OFFSET3P + 1);

    pIMiRNA3pSeq[0] = 'N';
    for (unsigned i = 0; i < length(pMiRNASeq); ++i) {
        pIMiRNASeq[i] = pMiRNASeq[i];
        if (i != 0 && i < SEED_REGION_LEN) {
            pIMiRNASeedSeq[idxseed] = pMiRNASeq[i];
            ++idxseed;
        }
        if (i >= SEED_REGION_LEN && i < length(pMiRNASeq) - OFFSET3P) {
            pIMiRNA3pSeq[idx3p] = pMiRNASeq[i];
            ++idx3p;
        }
    }
}

void MR3AlignScores::create_input_mrna_seq(
        mikan::TRNATYPE const &pMiRNASeq,
        mikan::TRNATYPE const &pMRNASeq,
        int pStart,
        int pEnd,
        const CharString &pSeedType,
        mikan::TRNATYPE &pIMRNASeq,
        mikan::TRNATYPE &pIMRNASeedSeq,
        Rna5String &pIMRNA3pSeq,
        bool &pNoMRNA1) {
    unsigned idx = 0;
    unsigned seqLen = (unsigned) (pEnd - pStart);
    int seed_idx = 0;
    int idex3p = 0;
    int len3p;
    int mirLen3p;
    unsigned maxpPos3p;
    unsigned seedRegLen;

    if (pSeedType == "7mer_BT" || pSeedType == "8mer_BT") {
        seedRegLen = SEED_REGION_LEN + 1;
    } else {
        seedRegLen = SEED_REGION_LEN;
    }

    resize(pIMRNASeq, seqLen);
    resize(pIMRNASeedSeq, seedRegLen - 1);

    len3p = seqLen - seedRegLen - OFFSET3P;
    if (len3p < 0) {
        len3p = 0;
    }

    mirLen3p = (int) length(pMiRNASeq) - seedRegLen - OFFSET3P;
    if (mirLen3p > len3p) {
        maxpPos3p = seqLen;
        len3p = seqLen - seedRegLen;
    } else {
        maxpPos3p = seqLen - OFFSET3P;
    }
    resize(pIMRNA3pSeq, len3p + 1);

    pIMRNA3pSeq[0] = 'N';
    for (unsigned i = 0; i < seqLen; ++i) {
        idx = seqLen - i - 1;
        if (pStart + i < length(pMRNASeq)) {
            pIMRNASeq[idx] = pMRNASeq[pStart + i];
        } else {
            pIMRNASeq[idx] = 'A';
            pNoMRNA1 = true;
        }

        if (idx != 0 && idx < seedRegLen) {
            pIMRNASeedSeq[length(pIMRNASeedSeq) - seed_idx - 1] = pIMRNASeq[idx];
            ++seed_idx;
        }
        if (idx >= seedRegLen && idx < maxpPos3p) {
            pIMRNA3pSeq[length(pIMRNA3pSeq) - idex3p - 1] = pIMRNASeq[idx];
            ++idex3p;
        }
    }
}

void MR3AlignScores::print_input(mikan::TRNATYPE &pInputMiRNASeq, mikan::TRNATYPE &pInputMRNASeq) {
    std::cout << "miRNA seq:   " << length(pInputMiRNASeq) << "," << pInputMiRNASeq;
    std::cout << std::endl;
    std::cout << "mRNA seq:    " << length(pInputMRNASeq) << "," << pInputMRNASeq;
    std::cout << std::endl;
}

//
// MR3EnergyScores methods
//

void MR3EnergyScores::clear_scores() {
    clear(mEffectiveSites);
    clear(mEnScores);
}

int MR3EnergyScores::calc_scores(
        MR3SeedSites &pSeedSites,
        mikan::TRNATYPE const &pMiRNASeq,
        mikan::TRNASet const &) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    std::string inputSeq;
    float score;

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    resize(mEnScores, length(pSeedSites.mEffectiveSites));

    mVRws.preppare_fold((int) length(pSeedSites.mEffectiveSites));

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            mEffectiveSites[i] = false;
            continue;
        }

        create_input_seq(i, pMiRNASeq, inputSeq);
//        print_input(inputSeq);

        mVRws.calc_fold_energy(i, inputSeq);
//        mVRws.print_fold_ret_vals(i);

        score = (float) mVRws.get_fold_energy(i);
        if (mMaxEnergy < score) {
            mEffectiveSites[i] = false;
            pSeedSites.mEffectiveSites[i] = false;
        } else {
            mEnScores[i] = score;
            mEffectiveSites[i] = true;
        }

    }

    return 0;
}

void
MR3EnergyScores::create_input_seq(int pIdx, mikan::TRNATYPE const &pMiRNASeq, std::string &pInputMRNASeq) {
    seqan::CharString inputMRNA;

    mAlign.get_mrna_seq(pIdx, inputMRNA);

    pInputMRNASeq.resize(length(pMiRNASeq) + LINKER_LEN + length(inputMRNA));

    int idx = 0;
    for (unsigned i = 0; i < length(pMiRNASeq); ++i) {
        pInputMRNASeq[idx] = pMiRNASeq[length(pMiRNASeq) - i - 1];
        ++idx;
    }
    for (unsigned i = 0; i < LINKER_LEN; ++i) {
        pInputMRNASeq[idx] = 'X';
        ++idx;
    }
    for (unsigned i = 0; i < length(inputMRNA); ++i) {
        pInputMRNASeq[idx] = inputMRNA[i];
        ++idx;
    }
}

void MR3EnergyScores::print_input(std::string &pInputSeq) {
    std::cout << "input seq: " << pInputSeq << std::endl;
}

//
// MR3SiteScores methods
//
void MR3SiteScores::clear_scores() {
    clear(mEffectiveSites);
    mAlignScores.clear_scores();
    mEnergyScores.clear_scores();
    mAlign.clear_align();
}

void MR3SiteScores::init_rnafold() {
    mVRws.set_backtrack(false);
    mVRws.set_fold_constraint(false);
    mVRws.init_arrays(RNAFOLD_MAX_INPUTLEN);
}

int MR3SiteScores::calc_scores(
        MR3SeedSites &pSeedSites,
        mikan::TRNATYPE const &miRNASeq,
        mikan::TRNASet const &pMRNASeqs) {
    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));

    mAlign.resize_align((int) length(pSeedSites.mEffectiveSites));

    mAlignScores.calc_scores(pSeedSites, miRNASeq, pMRNASeqs);
    mEnergyScores.calc_scores(pSeedSites, miRNASeq, pMRNASeqs);

    for (unsigned i = 0; i < length(mEffectiveSites); ++i) {
        if (mAlignScores.mEffectiveSites[i] && mEnergyScores.mEffectiveSites[i]) {
            mEffectiveSites[i] = true;
        } else {
            mEffectiveSites[i] = false;
        }
    }

    return 0;
}

void MR3SiteScores::print_alignment(int pIdx) {
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
// MR3TotalScores methods
//
void MR3TotalScores::clear_scores() {
    clear(mMRNAPos);
    clear(mSiteNum);
    clear(mTotalEnScores);
    clear(mTotalAlignScores);
    clear(mLogMaxAlignScores);
    clear(mLogMaxEnScores);
}

int MR3TotalScores::calc_scores(
        MR3SeedSites &pSeedSites,
        MR3SiteScores &pSiteScores,
        const seqan::String<unsigned> &pSortedSites) {

    const String<unsigned> &siteMRNAPos = pSeedSites.get_mrna_pos();
    int prevPos = -1;
    int newIdx = -1;
    String<int> newIdices;
    unsigned posIdx;
    String<double> maxAlignScores;
    String<double> maxEnScores;
    String<int> maxAlignScoreIds;
    String<int> maxEnScoreIds;
    String<bool> maxScoreProcessed;
    double scoreAlign;
    double scoreEn;
    double exp_diff;

    resize(newIdices, length(siteMRNAPos));
    for (unsigned i = 0; i < length(pSortedSites); ++i) {
        posIdx = pSortedSites[i];

        if (!pSiteScores.mEffectiveSites[posIdx]) {
            continue;
        }

        if (prevPos != (int) siteMRNAPos[posIdx]) {
            ++newIdx;
        }
        newIdices[i] = newIdx;
        prevPos = (int) siteMRNAPos[posIdx];
    }

    resize(mTotalAlignScores, newIdx + 1, 0.0);
    resize(mTotalEnScores, newIdx + 1, 0.0);
    resize(mLogMaxAlignScores, newIdx + 1, 0.0);
    resize(mLogMaxEnScores, newIdx + 1, 0.0);
    resize(mMRNAPos, newIdx + 1);
    resize(mSiteNum, newIdx + 1, 0);
    resize(maxAlignScores, newIdx + 1, 0.0);
    resize(maxAlignScoreIds, newIdx + 1);
    resize(maxEnScores, newIdx + 1, 0.0);
    resize(maxEnScoreIds, newIdx + 1);
    resize(maxScoreProcessed, newIdx + 1, false);

    for (unsigned i = 0; i < length(pSortedSites); ++i) {
        posIdx = pSortedSites[i];

        if (!pSiteScores.mEffectiveSites[posIdx]) {
            continue;
        }

        if (!maxScoreProcessed[newIdices[i]]) {
            maxScoreProcessed[newIdices[i]] = true;
            mLogMaxAlignScores[newIdices[i]] = 1.0;
            mLogMaxEnScores[newIdices[i]] = 1.0;
            mSiteNum[newIdices[i]] = 1;
            mMRNAPos[newIdices[i]] = siteMRNAPos[posIdx];
        } else {
            mSiteNum[newIdices[i]] += 1;
        }

        scoreAlign = pSiteScores.get_align_score(posIdx);
        mTotalAlignScores[newIdices[i]] += scoreAlign;
        if (scoreAlign > maxAlignScores[newIdices[i]]) {
            maxAlignScores[newIdices[i]] = scoreAlign;
            maxAlignScoreIds[newIdices[i]] = i;
        }

        scoreEn = -1.0 * pSiteScores.get_energy_score(posIdx);
        mTotalEnScores[newIdices[i]] += -1.0 * scoreEn;
        if (scoreEn > maxEnScores[newIdices[i]]) {
            maxEnScores[newIdices[i]] = scoreEn;
            maxEnScoreIds[newIdices[i]] = i;
        }
    }

    for (unsigned i = 0; i < length(pSortedSites); ++i) {
        posIdx = pSortedSites[i];

        if (!pSiteScores.mEffectiveSites[posIdx]) {
            continue;
        }

        if (maxAlignScores[newIdices[i]] != (int) i) {
            scoreAlign = pSiteScores.get_align_score(posIdx);
            exp_diff = scoreAlign - maxAlignScores[newIdices[i]];
            if (exp_diff > MIN_EXP_DIFF) {
                mLogMaxAlignScores[newIdices[i]] += std::exp(exp_diff);
            }
        }

        if (maxEnScores[newIdices[i]] != (int) i) {
            scoreEn = pSiteScores.get_energy_score(posIdx);
            exp_diff = scoreEn - maxEnScores[newIdices[i]];
            if (exp_diff > MIN_EXP_DIFF) {
                mLogMaxEnScores[newIdices[i]] += std::exp(exp_diff);
            }
        }

    }

    for (unsigned i = 0; i < length(maxAlignScores); ++i) {
        mLogMaxAlignScores[i] = maxAlignScores[i] + std::log(mTotalAlignScores[i]);
        mLogMaxEnScores[i] = 1.0 * (maxEnScores[i] + std::log(mTotalEnScores[i]));
    }

    return 0;
}

} // namespace mr3as
