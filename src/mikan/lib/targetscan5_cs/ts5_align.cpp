#include <sstream>               // stringstream
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_align.hpp"         // TS5Alignment

using namespace seqan;

namespace ts5cs {

//
// TS5Alignment methods
//
void TS5Alignment::clear_alignments() {
    clear(mAlignBars);
    clear(mAlignMRNA);
    clear(mAlignMiRNA);
}

void TS5Alignment::resize_alignments(unsigned pSize) {
    resize(mAlignBars, pSize);
    resize(mAlignMRNA, pSize);
    resize(mAlignMiRNA, pSize);
}

int TS5Alignment::align_seed(
        unsigned pMRNAIdx,
        CharString const &pSeedType,
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNAStr const &pMRNASeq,
        unsigned pSitePos) {
    int startUTR = 0;
    int seqLen = 0;
    int startMiRNA = 0;
    CharString matchBars;

    if (pSeedType == "8mer" || pSeedType == "7mer-m8") {
        startUTR = pSitePos - 1;
        seqLen = 8;
        startMiRNA = 7;
        matchBars = "||||||| ";
    } else if (pSeedType == "7mer-A1") {
        startUTR = pSitePos;
        seqLen = 7;
        startMiRNA = 6;
        matchBars = "|||||| ";
    }

    for (int i = 0; i < seqLen; ++i) {
        if (startUTR + i < (int) length(pMRNASeq)) {
            appendValue(mAlignMRNA[pMRNAIdx], pMRNASeq[startUTR + i]);
        } else {
            appendValue(mAlignMRNA[pMRNAIdx], ' ');
        }
    }

    for (int i = 0; i < seqLen; ++i) {
        appendValue(mAlignMiRNA[pMRNAIdx], pMiRNASeq[startMiRNA - i]);

    }

    for (unsigned i = 0; i < length(matchBars); ++i) {
        appendValue(mAlignBars[pMRNAIdx], matchBars[i]);
    }

    return 0;
}

int TS5Alignment::align_3p_part(
        unsigned pMRNAIdx,
        const CharString &,
        const mikan::TRNAStr &pMiRNAThreePrime,
        mikan::TRNAStr &pMRNAThreePrime,
        String<int> &pMatchLen,
        String<int> &pMiRNAPos,
        String<int> &pMRNAPos,
        float pScore,
        unsigned pMatchedIdx) {
    unsigned maxLen, lenDiff;
    mikan::TRNAStr miRNAThreePrime = pMiRNAThreePrime;

    maxLen = std::max(length(miRNAThreePrime), length(pMRNAThreePrime));
    if (pScore < 3) {
        return align_no_3p_part(pMRNAIdx, pMRNAThreePrime, miRNAThreePrime, maxLen);
    }

    complement(miRNAThreePrime);
    lenDiff = (unsigned) abs(pMRNAPos[pMatchedIdx] - pMiRNAPos[pMatchedIdx]);

    resize(mAlignBars[pMRNAIdx], maxLen + lenDiff, ' ');
    resize(mAlignMRNA[pMRNAIdx], maxLen + lenDiff, ' ');
    resize(mAlignMiRNA[pMRNAIdx], maxLen + lenDiff, ' ');

    set_align_bars(pMRNAIdx, pMatchLen, pMRNAPos, pMiRNAPos, pMatchedIdx);
    set_alignment(pMRNAIdx, pMRNAThreePrime, pMiRNAPos, pMRNAPos, pMatchedIdx, mAlignMRNA);
    set_alignment(pMRNAIdx, miRNAThreePrime, pMRNAPos, pMiRNAPos, pMatchedIdx, mAlignMiRNA);
    reverse(mAlignBars[pMRNAIdx]);
    reverse(mAlignMRNA[pMRNAIdx]);
    reverse(mAlignMiRNA[pMRNAIdx]);
    set_mismatch(pMRNAIdx, mAlignMRNA);
    set_mismatch(pMRNAIdx, mAlignMiRNA);

    trim_mismatch(pMRNAIdx);

    return 0;
}

int TS5Alignment::align_no_3p_part(
        unsigned pMRNAIdx,
        mikan::TRNAStr &pMRNAThreePrime,
        mikan::TRNAStr &pMiRNAThreePrime,
        unsigned pAlignLen) {
    mikan::TRNAStr miRNAThreePrime = pMiRNAThreePrime;
    complement(miRNAThreePrime);

    resize(mAlignBars[pMRNAIdx], pAlignLen, ' ');
    resize(mAlignMRNA[pMRNAIdx], pAlignLen, ' ');
    resize(mAlignMiRNA[pMRNAIdx], pAlignLen, ' ');

    for (unsigned i = 0; i < length(pMRNAThreePrime); ++i) {
        mAlignMRNA[pMRNAIdx][pAlignLen - i - 1] = pMRNAThreePrime[i];
    }
    for (unsigned i = 0; i < length(pMiRNAThreePrime); ++i) {
        mAlignMiRNA[pMRNAIdx][pAlignLen - i - 1] = miRNAThreePrime[i];
    }

    return 0;
}

void TS5Alignment::set_align_bars(
        unsigned pMRNAIdx,
        String<int> &pMatchLen,
        String<int> &pMRNAPos,
        String<int> &pMiRNAPos,
        unsigned pMatchedIdx) {
    unsigned startPos;

    startPos = (unsigned) std::max(pMRNAPos[pMatchedIdx], pMiRNAPos[pMatchedIdx]);
    for (int i = 0; i < pMatchLen[pMatchedIdx]; ++i) {
        mAlignBars[pMRNAIdx][startPos + i] = '|';
    }
}

void TS5Alignment::set_alignment(
        unsigned pMRNAIdx,
        mikan::TRNAStr &pSeqThreePrime,
        String<int> &pSeqPos1,
        String<int> &pSeqPos2,
        unsigned pMatchedIdx,
        mikan::TCharSet &pAlignSeq) {
    unsigned startPos;

    if (pSeqPos1[pMatchedIdx] - pSeqPos2[pMatchedIdx] > 0) {
        startPos = (unsigned) pSeqPos1[pMatchedIdx] - pSeqPos2[pMatchedIdx];
    } else {
        startPos = 0;
    }
    for (unsigned i = 0; i < length(pSeqThreePrime); ++i) {
        pAlignSeq[pMRNAIdx][startPos + i] = pSeqThreePrime[i];
    }
}

void TS5Alignment::set_mismatch(unsigned pMRNAIdx, mikan::TCharSet &pAlignSeq) {
    for (unsigned i = length(pAlignSeq[pMRNAIdx]) - 1; i > 0; --i) {
        if (pAlignSeq[pMRNAIdx][i] == ' ') {
            pAlignSeq[pMRNAIdx][i] = '-';
        } else {
            break;
        }
    }
}

void TS5Alignment::trim_mismatch(unsigned pMRNAIdx) {
    unsigned minLen, mmCount;
    unsigned alignMNALen = length(mAlignMRNA[pMRNAIdx]);
    unsigned alignMiRNALen = length(mAlignMiRNA[pMRNAIdx]);

    minLen = std::min(alignMNALen, alignMiRNALen);

    mmCount = 0;
    for (unsigned i = 0; i < minLen; ++i) {
        if (mAlignMRNA[pMRNAIdx][alignMNALen - 1 - i] == '-'
            && mAlignMiRNA[pMRNAIdx][alignMiRNALen - 1 - i] == '-') {
            ++mmCount;
        } else {
            break;
        }
    }

    if (mmCount > 0) {
        erase(mAlignMRNA[pMRNAIdx], alignMNALen - mmCount, alignMNALen);
        erase(mAlignMiRNA[pMRNAIdx], alignMiRNALen - mmCount, alignMiRNALen);
        erase(mAlignBars[pMRNAIdx], length(mAlignBars[pMRNAIdx]) - mmCount, length(mAlignBars[pMRNAIdx]));
    }
}

void TS5Alignment::write_alignment(int pIdx) const {
    std::stringstream stream;

    stream << "mRNA   5' " << mAlignMRNA[pIdx] << " 3'";
    stream << std::endl;
    stream << "          " << mAlignBars[pIdx] << "   ";
    stream << std::endl;
    stream << "miRNA  3' " << mAlignMiRNA[pIdx] << " 5'";
    stream << std::endl;

    std::cout << stream.str();
}

} // namespace ts5cs
