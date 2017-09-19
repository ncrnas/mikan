#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mr3_site_score.hpp"    // MR3DDGScores

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
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites) {

    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    const String<int> &mmPos = pSeedSites.get_mismatched_pos();
    const mikan::TCharSet &seedTypes = pSeedSites.get_seed_types();

    mikan::TRNAStr iMiRNASeq;
    mikan::TRNAStr iMRNASeq;
    mikan::TRNAStr iMiRNASeedSeq;
    mikan::TRNAStr iMRNASeedSeq;
    mikan::TRNAStr iMiRNA3pSeq;
    mikan::TRNAStr iMRNA3pSeq;
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
            seqEnd = sitePos[i] + (mikan::SEEDLEN + 2);
            seqStart = seqEnd - (TARGET_SEQ_LEN + 1);
            mm = mmPos[i];
        } else {
            seqEnd = sitePos[i] + (mikan::SEEDLEN + 1);
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
        mAlign.align_3p(i, iMiRNA3pSeq, iMRNA3pSeq);

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
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNAStr &pIMiRNASeq,
        mikan::TRNAStr &pIMiRNASeedSeq,
        mikan::TRNAStr &pIMiRNA3pSeq) {
    int idxseed = 0;
    int idx3p = 0;

    resize(pIMiRNASeq, length(pMiRNASeq));
    resize(pIMiRNASeedSeq, SEED_REGION_LEN - 1);
    resize(pIMiRNA3pSeq, length(pMiRNASeq) - SEED_REGION_LEN - OFFSET3P);

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
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNAStr const &pMRNASeq,
        int pStart,
        int pEnd,
        const mikan::TCharStr &pSeedType,
        mikan::TRNAStr &pIMRNASeq,
        mikan::TRNAStr &pIMRNASeedSeq,
        mikan::TRNAStr &pIMRNA3pSeq,
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
    resize(pIMRNA3pSeq, len3p);

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

void MR3AlignScores::print_input(mikan::TRNAStr &pInputMiRNASeq, mikan::TRNAStr &pInputMRNASeq) {
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
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &,
        mikan::MKSeedSites &pSeedSites) {

    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    std::string inputSeq;
    float score;

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    resize(mEnScores, length(pSeedSites.mEffectiveSites));

    mVRws.prepare_fold((int) length(pSeedSites.mEffectiveSites));

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

void MR3EnergyScores::create_input_seq(int pIdx, mikan::TRNAStr const &pMiRNASeq, std::string &pInputMRNASeq) {
    mikan::TCharStr inputMRNA;

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
void MR3SiteScores::init_from_args() {
    set_min_align_score(mOpts.mMinAlignScore);
    set_max_energy(mOpts.mMaxEnergy);

    mVRws.set_backtrack(false);
    mVRws.set_fold_constraint(false);
    mVRws.init_arrays(RNAFOLD_MAX_INPUTLEN);
}

void MR3SiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    mAlignScores.clear_scores();
    mEnergyScores.clear_scores();
    mAlign.clear_align();
}

int MR3SiteScores::calc_scores(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &) {

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));

    mAlign.resize_align((int) length(pSeedSites.mEffectiveSites));

    mAlignScores.calc_scores(pMiRNASeq, pMRNASeqs, pSeedSites);
    mEnergyScores.calc_scores(pMiRNASeq, pMRNASeqs, pSeedSites);

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

} // namespace mr3as
