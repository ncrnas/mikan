#include "mr3_align.hpp"    // MR3Align
#include "mk_typedef.hpp"   // TRNATYPE

using namespace seqan;
using namespace mikan;

namespace mr3as {

//
// MR3Align methods
//
template<class TRNAString>
void MR3Align<TRNAString>::clear_align() {
    clear(mEffectiveSites);
    clear(mAlignMRNA);
    clear(mAlignBars);
    clear(mAlignMiRNA);
    clear(mAlignSeedScores);
    clear(mAlign3PScores);
    clear(mAlignScores);
    clear(mAlignSeedMiRNA);
    clear(mAlignSeedMRNA);
    clear(mAlign3pMRNA);
    clear(mAlign3pMiRNA);
    clear(mGapCount3pMiRNA);
    clear(mGapCount3pMRNA);
}

template<class TRNAString>
void MR3Align<TRNAString>::resize_align(unsigned pSize) {
    resize(mEffectiveSites, pSize, false);
    resize(mAlignMRNA, pSize);
    resize(mAlignBars, pSize);
    resize(mAlignMiRNA, pSize);
    resize(mAlignSeedScores, pSize);
    resize(mAlign3PScores, pSize);
    resize(mAlignScores, pSize);
    resize(mAlignSeedMiRNA, pSize);
    resize(mAlignSeedMRNA, pSize);
    resize(mAlign3pMRNA, pSize);
    resize(mAlign3pMiRNA, pSize);
    resize(mGapCount3pMiRNA, pSize);
    resize(mGapCount3pMRNA, pSize);
}

template<class TRNAString>
void MR3Align<TRNAString>::align_seed(
        int pIdx,
        TRNAString &pIMiRNASeedSeq,
        TRNAString &pIMRNASeedSeq,
        int pMMpos) {
    int score = 0;
    unsigned idx2 = 0;

    resize(mAlignSeedMiRNA[pIdx], length(pIMRNASeedSeq));
    resize(mAlignSeedMRNA[pIdx], length(pIMRNASeedSeq));

    for (unsigned i = 0; i < length(pIMRNASeedSeq); ++i) {
        if (pMMpos == -1) {
            idx2 = i;
        } else if (i == length(pIMiRNASeedSeq) - pMMpos - 1) {
            mAlignSeedMiRNA[pIdx][i] = '-';
            mAlignSeedMRNA[pIdx][i] = pIMRNASeedSeq[i];
            score -= 36;
            continue;
        }

        mAlignSeedMiRNA[pIdx][i] = pIMiRNASeedSeq[idx2];
        mAlignSeedMRNA[pIdx][i] = pIMRNASeedSeq[i];
        if ((pIMiRNASeedSeq[idx2] == 'A' && pIMRNASeedSeq[i] == 'U')
            || (pIMiRNASeedSeq[idx2] == 'U' && pIMRNASeedSeq[i] == 'A')
            || (pIMiRNASeedSeq[idx2] == 'G' && pIMRNASeedSeq[i] == 'C')
            || (pIMiRNASeedSeq[idx2] == 'C' && pIMRNASeedSeq[i] == 'G')) {
            score += 20;
        } else if ((pIMiRNASeedSeq[idx2] == 'G' && pIMRNASeedSeq[i] == 'U')
                   || (pIMiRNASeedSeq[idx2] == 'U' && pIMRNASeedSeq[i] == 'G')) {
            score += 4;
        } else if (i != 0) {
            score -= 12;
        }

        if (pMMpos != -1) {
            ++idx2;
        }
    }

    mAlignSeedScores[pIdx] = score;

}

template<class TRNAString>
void MR3Align<TRNAString>::init_3p_align(int pIdx) {
    mAlign3PScores[pIdx] = 0;
    resize(mAlign3pMiRNA[pIdx], 0);
    resize(mAlign3pMRNA[pIdx], 0);
    mGapCount3pMiRNA[pIdx] = 0;
    mGapCount3pMRNA[pIdx] = 0;
}

template<class TRNAString>
void MR3Align<TRNAString>::align_3p(int pIdx, seqan::Rna5String &pIMiRNA3pSeq, seqan::Rna5String &pIMRNA3pSeq) {
    int score;

    clearClipping(mAign3P);
    clearGaps(mAign3P);
    resize(rows(mAign3P), 2);
    assignSource(row(mAign3P, 0), pIMiRNA3pSeq);
    assignSource(row(mAign3P, 1), pIMRNA3pSeq);

    score = localAlignment(mAign3P, mScoreMatrix3P);
    score -= 140;

//    std::cout << score << std::endl;
//    std::cout << mAign3P << std::endl;

    if (score > 0) {
        mAlign3PScores[pIdx] = score;

        TGap &alignMiRNA = (TGap &) row(mAign3P, 0);
        TGap &alignMRNA = (TGap &) row(mAign3P, 1);

        resize(mAlign3pMiRNA[pIdx], length(alignMiRNA) - 1);
        resize(mAlign3pMRNA[pIdx], length(alignMRNA) - 1);

        for (unsigned i = 0; i < length(alignMiRNA) - 1; ++i) {
            mAlign3pMiRNA[pIdx][i] = alignMiRNA[i + 1];
            mAlign3pMRNA[pIdx][i] = alignMRNA[i + 1];
            if (alignMiRNA[i + 1] == '-') {
                ++mGapCount3pMiRNA[pIdx];
            }
            if (alignMRNA[i + 1] == '-') {
                ++mGapCount3pMRNA[pIdx];
            }
        }
    }
}

template<class TRNAString>
void MR3Align<TRNAString>::combine_alignments(
        int pIdx,
        TRNAString const &pMiRNASeq,
        TRNAString const &pMRNASeq,
        bool noA1) {
    int maxlen = (int) length(pMiRNASeq) + mGapCount3pMiRNA[pIdx];
    int idx2;
    int idx2_orig;

    resize(mAlignMiRNA[pIdx], maxlen);
    resize(mAlignMRNA[pIdx], maxlen);
    resize(mAlignBars[pIdx], maxlen);

    mAlignMiRNA[pIdx][0] = pMiRNASeq[0];
    if (noA1) {
        mAlignMRNA[pIdx][0] = '-';
    } else {
        mAlignMRNA[pIdx][0] = pMRNASeq[0];
    }

    for (unsigned i = 0; i < length(mAlignSeedMiRNA[pIdx]); ++i) {
        idx2 = i + 1;
        mAlignMiRNA[pIdx][idx2] = mAlignSeedMiRNA[pIdx][i];
        mAlignMRNA[pIdx][idx2] = mAlignSeedMRNA[pIdx][i];
    }
    for (unsigned i = 0; i < length(mAlign3pMiRNA[pIdx]); ++i) {
        idx2 = i + (int) length(mAlignSeedMiRNA[pIdx]) + 1;
        mAlignMiRNA[pIdx][idx2] = mAlign3pMiRNA[pIdx][i];

        if (i < length(mAlignMRNA[pIdx]) - 1) {
            mAlignMRNA[pIdx][idx2] = mAlign3pMRNA[pIdx][i];
        } else {
            mAlignMRNA[pIdx][idx2] = '-';
            ++mGapCount3pMRNA[pIdx];
        }
    }

    idx2 = (int) length(mAlignSeedMiRNA[pIdx]) + (int) length(mAlign3pMiRNA[pIdx]) + 1;
    idx2_orig = idx2 - mGapCount3pMiRNA[pIdx];
    while (idx2 < maxlen) {
        if (idx2_orig < (int) length(pMiRNASeq)) {
            mAlignMiRNA[pIdx][idx2] = pMiRNASeq[idx2_orig];
        } else {
            mAlignMiRNA[pIdx][idx2] = '-';
            ++mGapCount3pMiRNA[pIdx];
        }
        ++idx2;
        ++idx2_orig;
    }

    idx2 = (int) length(mAlignSeedMiRNA[pIdx]) + (int) length(mAlign3pMiRNA[pIdx]) + 1;
    idx2_orig = idx2 - mGapCount3pMRNA[pIdx];
    while (idx2 < maxlen) {
        if (idx2_orig < (int) length(pMRNASeq)) {
            mAlignMRNA[pIdx][idx2] = pMRNASeq[idx2_orig];
        } else {
            mAlignMRNA[pIdx][idx2] = '-';
            ++mGapCount3pMRNA[pIdx];
        }
        ++idx2;
        ++idx2_orig;
    }

    set_align_bars(pIdx);

    reverse(mAlignMiRNA[pIdx]);
    reverse(mAlignMRNA[pIdx]);
    reverse(mAlignBars[pIdx]);

//    std::cout <<  mAlignMiRNA[pIdx] << std::endl;
//    std::cout <<  mAlignBars[pIdx] << std::endl;
//    std::cout << mAlignMRNA[pIdx] << std::endl;
}

template<class TRNAString>
void MR3Align<TRNAString>::set_align_bars(int pIdx) {
    char miRNA, mRNA;

    for (unsigned i = 0; i < length(mAlignBars[pIdx]); ++i) {
        miRNA = mAlignMiRNA[pIdx][i];
        mRNA = mAlignMRNA[pIdx][i];

        if (miRNA == 'A') {
            if (mRNA == 'U') {
                mAlignBars[pIdx][i] = '|';
            } else {
                mAlignBars[pIdx][i] = ' ';
            }
        } else if (miRNA == 'C') {
            if (mRNA == 'G') {
                mAlignBars[pIdx][i] = '|';
            } else {
                mAlignBars[pIdx][i] = ' ';
            }
        } else if (miRNA == 'G') {
            if (mRNA == 'C') {
                mAlignBars[pIdx][i] = '|';
            } else if (mRNA == 'U') {
                mAlignBars[pIdx][i] = ':';
            } else {
                mAlignBars[pIdx][i] = ' ';
            }
        } else if (miRNA == 'U') {
            if (mRNA == 'A') {
                mAlignBars[pIdx][i] = '|';
            } else if (mRNA == 'G') {
                mAlignBars[pIdx][i] = ':';
            } else {
                mAlignBars[pIdx][i] = ' ';
            }
        } else {
            mAlignBars[pIdx][i] = ' ';
        }
    }

}

template<class TRNAString>
void MR3Align<TRNAString>::get_mrna_seq(int pIdx, seqan::CharString &pStrMRNA) {
    int idx = 0;

    resize(pStrMRNA, length(mAlignMRNA[pIdx]) - mGapCount3pMRNA[pIdx]);
    for (unsigned i = 0; i < length(mAlignMRNA[pIdx]); ++i) {
        if (mAlignMRNA[pIdx][i] != '-' || i == length(mAlignMRNA[pIdx]) - 1) {
            pStrMRNA[idx] = mAlignMRNA[pIdx][i];
            ++idx;
        }
    }

    resize(pStrMRNA, idx);
    reverse(pStrMRNA);
}

// Explicit template instantiation
template
class MR3Align<TRNATYPE>;

} // namespace mr3as
