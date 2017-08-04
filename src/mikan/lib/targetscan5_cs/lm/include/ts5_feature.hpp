#ifndef TS5_FEATURE_HPP_
#define TS5_FEATURE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_align.hpp"         // TS5Alignment
#include "ts5_seed_site.hpp"     // TS5SeedSites

namespace ts5cs {

//
// Seed site position feature
//
class TS5FeatSitePos {
public:
    // Define methods
    TS5FeatSitePos() : mMaxLen(1500) {}

    int get_val(int i) { return mSitePos[i]; }

    // Method prototype
    int add_features(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     mikan::MKSeedSites &pSeedSites);

    void clear_features();

private:
    // Define variables
    seqan::String<int> mSitePos;
    const int mMaxLen;

};

//
// AU-rich feature
//
class TS5FeatAURich {
public:
    // Define methods
    TS5FeatAURich() {}

    float get_val(int i) { return mAURich[i]; }

    // Method prototype
    int add_features(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     mikan::MKSeedSites &pSeedSites);

    void clear_features();

private:
    seqan::String<float> mAURich;

private:
    void getUpDownStreamPos(seqan::CharString pSeedType, int pStartPos, int &pStartU, int &pStartD);

    void calcPosScores(const seqan::CharString &pSeedType, seqan::CharString &pUpOrDown,
                       const mikan::TRNAStr &pAURichRNA, int pStart, int pEnd, float &pTotalScore, float &pMaxScore);
};

//
//  Additional 3' paring feature
//
class TS5FeatThreePrimePair {
public:
    // Define methods
    TS5FeatThreePrimePair() : mIdxBestScore(0) {}

    float get_val(int i) { return mThreePrimePair[i]; }

    const TS5Alignment &get_alignment() { return mAlign; }

    // Method prototype
    int add_features(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     mikan::MKSeedSites &pSeedSites);

    void clear_features();

private:
    seqan::String<float> mThreePrimePair;
    unsigned mIdxBestScore;
    TS5Alignment mAlign;

private:
    void getUpDownStreamPos(seqan::CharString pSeedType, int pStartPos, int &pStart, int &pEnd);

    void getMRNASeq(const seqan::CharString &pSeedType, const mikan::TRNAStr &pMRNASeq, int pSitePos,
                    mikan::TRNAStr &pMRNAThreePrime);

    void getMiRNASeq(const seqan::CharString &pSeedType, const mikan::TRNAStr &pMiRNASeq,
                     mikan::TRNAStr &pMiRNAThreePrime);

    float findBestMatch(unsigned pPosIdx, mikan::TSitePosSet const &pMRNAPos, mikan::TSitePosSet const &pSitePos,
                        const seqan::CharString &pSeedType, const mikan::TRNAStr &pMRNASeq,
                        const mikan::TRNAStr &pMiRNASeq);

    void connectMatchedSeq(seqan::String<int> &pMatchLen, seqan::String<int> &pMiRNAPos,
                           seqan::String<int> &pMRNAPos);

    float calcScore(const seqan::CharString &pSeedType, seqan::String<int> &pMatchLen,
                    seqan::String<int> &pMiRNAPos, seqan::String<int> &pMRNAPos);

};

//
// Store all raw feature values
//
class TS5RawFeatures {
public:
    // Define methods
    TS5RawFeatures() {}

    int get_site_pos(int i) { return mSitePos.get_val(i); }

    float get_au_rich(int i) { return mAURich.get_val(i); }

    float get_three_prime_pair(int i) { return mThreePrimePair.get_val(i); }

    const TS5Alignment &get_alignment() { return mThreePrimePair.get_alignment(); }

    // Method prototypes
    int add_features(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     mikan::MKSeedSites &pSeedSites);

    void clear_features();

private:
    TS5FeatSitePos mSitePos;
    TS5FeatAURich mAURich;
    TS5FeatThreePrimePair mThreePrimePair;
};

} // namespace ts5cs

#endif /* TS5_FEATURE_HPP_ */

