#include <seqan/seq_io.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"        // MKSeedSeqs

using namespace seqan;

namespace mikan {

//
// MKSeedSeqs methods
//
void MKSeedSeqs::init_from_args() {
    resize(mSeedTypeDef, 1);
    mSeedTypeDef[0] = "";
}

void MKSeedSeqs::set_flags() {
    mSingleGU = true;
    mMultiGU = true;
    mMisMatch = true;
    mGUMisMatch = true;
    mBT = true;
    mBTM8 = true;
    mBM = true;
    mLP = true;
    mOther = true;
    mAddInReverse = false;
}

void MKSeedSeqs::clear_seeds() {
    nNumNewSeq = 0;

    clear(mMiRNASeq);

    clear(mEffectiveSeeds);
    clear(mSeedSeqs);
    clear(mSeedTypes);
    clear(mMisMatchPos);
}

int MKSeedSeqs::create_seed_seqs(mikan::TRNAStr const &pSeq) {
    if (length(pSeq) == 0) {
        return 1;
    }
    mMiRNASeq = pSeq;

    mikan::TRNAStr seedSeq;
    CharString seedType;

    int retVal;

    resize(seedSeq, 6);
    for (unsigned i = 0; i < length(seedSeq); ++i) {
        seedSeq[i] = mMiRNASeq[i + 1];
    }
    reverseComplement(seedSeq);

    appendValue(mSeedSeqs, seedSeq);
    appendValue(mSeedTypes, mNMerLab);
    appendValue(mMisMatchPos, 0);

    if (mSingleGU) {
        retVal = create_single_guwobble_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mMultiGU) {
        retVal = create_multi_guwobble_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mMisMatch) {
        retVal = create_mismatch_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mGUMisMatch) {
        retVal = create_gu_mismatch_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mBT) {
        retVal = create_bt_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mBTM8) {
        retVal = create_bt_m8_seed_seqs();
        if (retVal != 0) {
            return 1;
        }
    }

    if (mBM) {
        retVal = create_bm_seed_seqs();
        if (retVal != 0) {
            return 1;
        }
    }

    if (mLP) {
        mMMLab = "LP";
        retVal = create_mismatch_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mOther) {
        retVal = create_other_seed_seqs(seedSeq);
        if (retVal != 0) {
            return 1;
        }
    }

    resize(mEffectiveSeeds, length(mSeedSeqs), true);

    if (mFilterRedundant) {
        retVal = check_redundant_seeds();
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;
}

int MKSeedSeqs::create_single_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq) {
    mikan::TRNAStr seedGUSeq;
    nNumNewSeq = 0;
    unsigned mmPos;

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        if (mTSSVMMismatch) {
            mmPos = length(pSeedSeq) - i;
        } else {
            mmPos = i;
        }

        if (pSeedSeq[i] == 'C') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            mTmpSeedSeqs[nNumNewSeq] = seedGUSeq;
            mTmpSeedTypes[nNumNewSeq] = mGUTLab;
            mTmpMisMatchPos[nNumNewSeq] = mmPos;
            nNumNewSeq++;
        } else if (pSeedSeq[i] == 'A') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            mTmpSeedSeqs[nNumNewSeq] = seedGUSeq;
            mTmpSeedTypes[nNumNewSeq] = mGUMLab;
            mTmpMisMatchPos[nNumNewSeq] = mmPos;
            nNumNewSeq++;
        }
    }

    int retVal = add_seed_seqs();

    return retVal;
}

int MKSeedSeqs::create_multi_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq) {
    mikan::TRNAStr seedGUSeq;
    nNumNewSeq = 0;
    unsigned mm;

    unsigned baseLen = length(mSeedSeqs);
    unsigned vecSize;

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        vecSize = baseLen + nNumNewSeq;
        if (pSeedSeq[i] == 'C') {
            for (unsigned j = 1; j < vecSize; ++j) {
                if (j < baseLen) {
                    seedGUSeq = get_seed_seq(j);
                    mm = get_mismatched_pos(j);
                } else {
                    seedGUSeq = mTmpSeedSeqs[j - baseLen];
                    mm = mTmpMisMatchPos[j - baseLen];
                }
                if (seedGUSeq[i] != 'U' && mm > i) {
                    seedGUSeq[i] = 'U';
                    mTmpSeedSeqs[nNumNewSeq] = seedGUSeq;
                    mTmpSeedTypes[nNumNewSeq] = mMultiGUTLab;
                    mTmpMisMatchPos[nNumNewSeq] = mm;
                    nNumNewSeq++;
                }
            }
        } else if (pSeedSeq[i] == 'A') {
            for (unsigned j = 1; j < vecSize; ++j) {
                if (j < baseLen) {
                    seedGUSeq = get_seed_seq(j);
                    mm = get_mismatched_pos(j);
                } else {
                    seedGUSeq = mTmpSeedSeqs[j - baseLen];
                    mm = mTmpMisMatchPos[j - baseLen];
                }
                if (seedGUSeq[i] != 'G' && mm > i) {
                    seedGUSeq[i] = 'G';
                    mTmpSeedSeqs[nNumNewSeq] = seedGUSeq;
                    mTmpSeedTypes[nNumNewSeq] = mMultiGUMLab;
                    mTmpMisMatchPos[nNumNewSeq] = mm;
                    nNumNewSeq++;
                }
            }
        }
    }

    int retVal = add_seed_seqs();

    return retVal;
}

int MKSeedSeqs::create_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq, bool pIsMMGU, int pGUPos) {
    mikan::TRNAStr seedLPSeq;
    CharString seedType;
    nNumNewSeq = 0;
    char ch1 = 0;
    char ch2 = 0;
    char ch3 = 0;
    unsigned mmPos;

    if (pIsMMGU) {
        seedType = mMMGULab;
    } else {
        seedType = mMMLab;
    }

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        if (pIsMMGU && (pGUPos == static_cast<int>(i))) {
            continue;
        }
        if (pSeedSeq[i] == 'A') {
            ch1 = 'C';
            ch2 = 'x';
            ch3 = 'U';
        } else if (pSeedSeq[i] == 'C') {
            ch1 = 'A';
            ch2 = 'x';
            ch3 = 'G';
        } else if (pSeedSeq[i] == 'G') {
            ch1 = 'A';
            ch2 = 'U';
            ch3 = 'C';
        } else if (pSeedSeq[i] == 'U') {
            ch1 = 'C';
            ch2 = 'G';
            ch3 = 'A';
        }

        if (mTSSVMMismatch) {
            mmPos = length(pSeedSeq) - i;
        } else {
            mmPos = i;
        }

        seedLPSeq = pSeedSeq;
        seedLPSeq[i] = ch1;
        mTmpSeedSeqs[nNumNewSeq] = seedLPSeq;
        mTmpSeedTypes[nNumNewSeq] = seedType;
        mTmpMisMatchPos[nNumNewSeq] = mmPos;
        nNumNewSeq++;
        if (ch2 != 'x') {
            seedLPSeq = pSeedSeq;
            seedLPSeq[i] = ch2;
            mTmpSeedSeqs[nNumNewSeq] = seedLPSeq;
            mTmpSeedTypes[nNumNewSeq] = seedType;
            mTmpMisMatchPos[nNumNewSeq] = mmPos;
            nNumNewSeq++;
        }
        seedLPSeq = pSeedSeq;
        seedLPSeq[i] = ch3;
        mTmpSeedSeqs[nNumNewSeq] = seedLPSeq;
        mTmpSeedTypes[nNumNewSeq] = seedType;
        mTmpMisMatchPos[nNumNewSeq] = mmPos;
        nNumNewSeq++;
    }

    int retVal = add_seed_seqs();

    return retVal;
}

int MKSeedSeqs::create_gu_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq) {
    mikan::TRNAStr seedGUSeq;
    int retVal;

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        retVal = 0;
        if (pSeedSeq[i] == 'C') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            retVal = create_mismatch_seed_seqs(seedGUSeq, true, i);
        } else if (pSeedSeq[i] == 'A') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            retVal = create_mismatch_seed_seqs(seedGUSeq, true, i);
        }
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;
}

int MKSeedSeqs::create_bt_seed_seqs(mikan::TRNAStr &pSeedSeq) {
    mikan::TRNAStr seedBTSeq;
    nNumNewSeq = 0;

    resize(seedBTSeq, 6);
    for (unsigned i = 0; i < length(seedBTSeq); ++i) {
        for (unsigned j = 0; j < length(mRNAChar); ++j) {
            int l = 0;
            for (unsigned k = 0; k < length(seedBTSeq); ++k) {
                if (i == k) {
                    seedBTSeq[k] = mRNAChar[j];
                    mTmpMisMatchPos[nNumNewSeq] = i;
                } else {
                    seedBTSeq[k] = pSeedSeq[l];
                    ++l;
                }
            }
            mTmpSeedSeqs[nNumNewSeq] = seedBTSeq;
            mTmpSeedTypes[nNumNewSeq] = "BT";
            nNumNewSeq++;
        }
    }

    int retVal = add_seed_seqs();

    return retVal;
}

int MKSeedSeqs::create_bt_m8_seed_seqs() {
    mikan::TRNAStr seedBTSeq;
    nNumNewSeq = 0;

    resize(seedBTSeq, 6);
    for (unsigned i = 0; i < length(seedBTSeq); ++i) {
        for (unsigned j = 0; j < length(mRNAChar); ++j) {
            int l = 0;
            for (unsigned k = 0; k < length(seedBTSeq); ++k) {
                if (i == k) {
                    seedBTSeq[k] = mRNAChar[j];
                } else {
                    seedBTSeq[k] = mMiRNASeq[l + 2];
                    ++l;
                }
            }
            reverseComplement(seedBTSeq);
            mTmpSeedSeqs[nNumNewSeq] = seedBTSeq;
            mTmpSeedTypes[nNumNewSeq] = "BT";
            mTmpMisMatchPos[nNumNewSeq] = i + 2;
            nNumNewSeq++;
        }
    }

    int retVal = add_seed_seqs();

    return retVal;
}

int MKSeedSeqs::create_bm_seed_seqs() {
    mikan::TRNAStr seedBMSeq;
    nNumNewSeq = 0;

    resize(seedBMSeq, 6);
    for (unsigned i = 0; i < length(seedBMSeq); ++i) {
        int k = 0;
        for (unsigned j = 0; j < length(seedBMSeq) + 1; ++j) {
            if (i == j) {
                continue;
            }
            seedBMSeq[k] = mMiRNASeq[j + 1];
            ++k;
        }
        reverseComplement(seedBMSeq);
        mTmpSeedSeqs[nNumNewSeq] = seedBMSeq;
        mTmpSeedTypes[nNumNewSeq] = "BM";
        mTmpMisMatchPos[nNumNewSeq] = i + 1;
        nNumNewSeq++;
    }

    int retVal = add_seed_seqs();

    return retVal;
}

int MKSeedSeqs::create_other_seed_seqs(mikan::TRNAStr &) {
    return 0;
}

int MKSeedSeqs::check_redundant_seeds() {
    mikan::TRNAStr seedSeq;
    mikan::TIndexQGram RNAIdx(mSeedSeqs);
    mikan::TFinder finder(RNAIdx);

    for (unsigned i = 0; i < length(mSeedSeqs); ++i) {
        if (!mEffectiveSeeds[i]) {
            continue;
        }
        seedSeq = mSeedSeqs[i];
        while (find(finder, seedSeq)) {
            unsigned pos = position(finder).i1;
            if (i != pos
                && ((mSeedTypes[i] != "MM" && mSeedTypes[i] != "MMGU" && mSeedTypes[i] != "BT")
                    || (mSeedTypes[i] == "MM" && mSeedTypes[pos] == "MM")
                    || (mSeedTypes[i] == "MM" && mSeedTypes[pos] == "MMGU")
                    || (mSeedTypes[i] == "MMGU" && mSeedTypes[pos] == "MMGU")
                    || (mSeedTypes[i] == "BT" && mSeedTypes[pos] == "BT"))) {
                mEffectiveSeeds[pos] = false;
            }
        }
        goBegin(finder);
        clear(finder);
    }

    return 0;
}

void MKSeedSeqs::init_temp(unsigned pVecSize) {
    nNumNewSeq = 0;
    resize(mTmpSeedSeqs, pVecSize);
    resize(mTmpSeedTypes, pVecSize);
    resize(mTmpMisMatchPos, pVecSize);
}

int MKSeedSeqs::add_seed_seqs() {

    unsigned idxTemp;
    unsigned baseSize = length(mSeedSeqs);

    resize(mSeedSeqs, baseSize + nNumNewSeq);
    resize(mSeedTypes, baseSize + nNumNewSeq);
    resize(mMisMatchPos, baseSize + nNumNewSeq);

    for (unsigned i = 0; i < nNumNewSeq; ++i) {
        if (mAddInReverse) {
            idxTemp = nNumNewSeq - i - 1;
        } else {
            idxTemp = i;
        }
        mSeedSeqs[baseSize + i] = mTmpSeedSeqs[idxTemp];
        mSeedTypes[baseSize + i] = mTmpSeedTypes[idxTemp];
        mMisMatchPos[baseSize + i] = mTmpMisMatchPos[idxTemp];
    }

    return 0;
}

void MKSeedSeqs::print_all() {
    unsigned numSeqs = length(mSeedSeqs);
    mikan::TRNAStr seedRC;

    std::cout << "The number of seeds: " << numSeqs << std::endl;

    for (unsigned i = 0; i < numSeqs; i++) {
        seedRC = mSeedSeqs[i];
        reverseComplement(seedRC);
        std::cout << i << ", ";
        std::cout << mSeedSeqs[i] << ", " << seedRC << ", ";
        std::cout << mSeedTypes[i] << ", ";
        std::cout << mMisMatchPos[i] << ", ";
        std::cout << mEffectiveSeeds[i] << std::endl;
    }
}

} // namespace mikan
