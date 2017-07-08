#ifndef TS5_FEATURE_HPP_
#define TS5_FEATURE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_align.hpp"         // TS5Alignment
#include "ts5_seed_site.hpp"     // TS5SeedSites

namespace ts5cs {

//
// Seed type feature
//
class TS5FeatSeedType {
public:
    // Define methods
    TS5FeatSeedType() {}

    seqan::CharString get_seed_type(unsigned idx) { return mSeedTypes[idx]; }

    seqan::CharString &get_val(int i) { return mSeedTypes[i]; }

    // Method prototype
    int add_features(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     seqan::String<bool> &pEffectiveSites, mikan::TSitePos const &pMRNAPos, mikan::TSitePos const &pSitePos);

    void clear_features();

private:
    seqan::StringSet<seqan::CharString> mSeedTypes;
};

//
// Seed site position feature
//
class TS5FeatSitePos {
public:
    // Constant values
    static const int MIN_DIST_TO_CDS = 15;

public:
    // Define methods
    TS5FeatSitePos() : mMaxLen(1500) {}

    int &get_val(int i) { return mSitePos[i]; }

    // Method prototype
    int add_features(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     seqan::String<bool> &pEffectiveSites, mikan::TSitePos const &pMRNAPos, mikan::TSitePos const &pSitePos,
                     TS5FeatSeedType &pSeedTypes);

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

    float &get_val(int i) { return mAURich[i]; }

    // Method prototype
    int add_features(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     seqan::String<bool> &pEffectiveSites, mikan::TSitePos const &pMRNAPos, mikan::TSitePos const &pSitePos,
                     TS5FeatSeedType &pSeedTypes);

    void clear_features();

private:
    seqan::String<float> mAURich;

private:
    void getUpDownStreamPos(seqan::CharString pSeedType, int pStartPos, int &pStartU, int &pStartD);

    void calcPosScores(const seqan::CharString &pSeedType, seqan::CharString &pUpOrDown,
                       const mikan::TRNATYPE &pAURichRNA, int pStart, int pEnd, float &pTotalScore, float &pMaxScore);
};

//
//  Additional 3' paring feature
//
class TS5FeatThreePrimePair {
public:
    // Define methods
    TS5FeatThreePrimePair() : mIdxBestScore(0) {}

    float &get_val(int i) { return mThreePrimePair[i]; }

    const TS5Alignment &get_alignment() { return mAlign; }

    // Method prototype
    int add_features(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     seqan::String<bool> &pEffectiveSites, mikan::TSitePos const &pMRNAPos, mikan::TSitePos const &pSitePos,
                     TS5FeatSeedType &pSeedTypes);

    void clear_features();

private:
    seqan::String<float> mThreePrimePair;
    unsigned mIdxBestScore;
    TS5Alignment mAlign;

private:
    void getUpDownStreamPos(seqan::CharString pSeedType, int pStartPos, int &pStart, int &pEnd);

    void getMRNASeq(const seqan::CharString &pSeedType, const mikan::TRNATYPE &pMRNASeq, int pSitePos,
                    mikan::TRNATYPE &pMRNAThreePrime);

    void getMiRNASeq(const seqan::CharString &pSeedType, const mikan::TRNATYPE &pMiRNASeq,
                     mikan::TRNATYPE &pMiRNAThreePrime);

    float findBestMatch(unsigned pPosIdx, mikan::TSitePos const &pMRNAPos, mikan::TSitePos const &pSitePos,
                        const seqan::CharString &pSeedType, const mikan::TRNATYPE &pMRNASeq,
                        const mikan::TRNATYPE &pMiRNASeq);

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
    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    TS5RawFeatures() {}

    seqan::CharString &get_seed_type(int i) { return mSeedTypes.get_val(i); }

    int &get_site_pos(int i) { return mSitePos.get_val(i); }

    float &get_au_rich(int i) { return mAURich.get_val(i); }

    float &get_three_prime_pair(int i) { return mThreePrimePair.get_val(i); }

    const TS5Alignment &get_alignment() { return mThreePrimePair.get_alignment(); }

    // Method prototypes
    int add_features(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                     TS5SeedSites const &pSeedSites);

    void clear_features();

private:
    TS5FeatSeedType mSeedTypes;
    TS5FeatSitePos mSitePos;
    TS5FeatAURich mAURich;
    TS5FeatThreePrimePair mThreePrimePair;
};

} // namespace ts5cs

#endif /* TS5_FEATURE_HPP_ */

