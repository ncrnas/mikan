#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_site_score.hpp"   // PITASiteScores

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
        mikan::TCharStr const &pSeedType,
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
    startPos = pSitePos + (mikan::SEEDLEN + 1) - (int) length(pMiRNASeq);
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

} // namespace ptddg
