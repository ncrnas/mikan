#ifndef TS5_FEATURE_HPP_
#define TS5_FEATURE_HPP_

#include <mikan/lib/targetscan5_cs/include/ts5_align.hpp>         // TS5Alignment
#include <mikan/lib/targetscan5_cs/include/ts5_seed_site.hpp>     // TS5SeedSites
#include <seqan/sequence.h>

namespace ts5cs{

//
// Seed type feature
//
template <class TRNAString>
class TS5FeatSeedType
{
public:
    // Define types
    typedef typename seqan::StringSet<TRNAString> TRNASet;
    typedef typename seqan::StringSet<seqan::CharString> TCharSet;
    typedef typename seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TS5FeatSeedType() {}
    seqan::CharString get_seed_type(unsigned idx){return mSeedTypes[idx];}
    seqan::CharString& get_val(int i){return mSeedTypes[i];}

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos, TSitePos const &pSitePos);
    void clear_features();

private:
    seqan::StringSet<seqan::CharString> mSeedTypes;
};

//
// Seed site position feature
//
template <class TRNAString>
class TS5FeatSitePos
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::String<unsigned> TSitePos;

    // Constant values
    static const int MIN_DIST_TO_CDS = 15;

public:
    // Define methods
    TS5FeatSitePos(): mMaxLen(1500) {}
    int& get_val(int i){return mSitePos[i];}

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos, TSitePos const &pSitePos,
            TS5FeatSeedType<TRNAString> &pSeedTypes);
    void clear_features();

private:
    // Define variables
    seqan::String<int> mSitePos;
    const int mMaxLen;

};

//
// AU-rich feature
//
template <class TRNAString>
class TS5FeatAURich
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TS5FeatAURich() {}
    float& get_val(int i){return mAURich[i];}

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos, TSitePos const &pSitePos,
            TS5FeatSeedType<TRNAString> &pSeedTypes);
    void clear_features();

private:
    seqan::String<float> mAURich;

private:
    void getUpDownStreamPos(seqan::CharString pSeedType, int pStartPos, int& pStartU, int& pStartD);
    void calcPosScores(const seqan::CharString& pSeedType, seqan::CharString& pUpOrDown,
            const TRNAString& pAURichRNA, int pStart, int pEnd, float& pTotalScore, float& pMaxScore);
};

//
//  Additional 3' paring feature
//
template <class TRNAString>
class TS5FeatThreePrimePair
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TS5FeatThreePrimePair(): mIdxBestScore(0) {}
    float& get_val(int i){return mThreePrimePair[i];}
    const TS5Alignment<TRNAString>& get_alignment() {return mAlign;}

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos, TSitePos const &pSitePos,
            TS5FeatSeedType<TRNAString> &pSeedTypes);
    void clear_features();

private:
    seqan::String<float> mThreePrimePair;
    unsigned mIdxBestScore;
    TS5Alignment<TRNAString> mAlign;

private:
    void getUpDownStreamPos(seqan::CharString pSeedType, int pStartPos, int& pStart, int& pEnd);
    void getMRNASeq(const seqan::CharString &pSeedType, const TRNAString &pMRNASeq, int pSitePos,
            TRNAString &pMRNAThreePrime);
    void getMiRNASeq(const seqan::CharString &pSeedType, const TRNAString &pMiRNASeq,
            TRNAString &pMiRNAThreePrime);
    float findBestMatch(unsigned pPosIdx, TSitePos const &pMRNAPos, TSitePos const &pSitePos,
            const seqan::CharString &pSeedType, const TRNAString &pMRNASeq, const TRNAString &pMiRNASeq);
    void connectMatchedSeq(seqan::String<int> &pMatchLen, seqan::String<int> &pMiRNAPos,
            seqan::String<int> &pMRNAPos);
    float calcScore(const seqan::CharString &pSeedType, seqan::String<int> &pMatchLen,
            seqan::String<int> &pMiRNAPos, seqan::String<int> &pMRNAPos);

};

//
// Store all raw feature values
//
template <class TRNAString>
class TS5RawFeatures
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    TS5RawFeatures() {}
    seqan::CharString& get_seed_type(int i){return mSeedTypes.get_val(i);}
    int& get_site_pos(int i){return mSitePos.get_val(i);}
    float& get_au_rich(int i){return mAURich.get_val(i);}
    float& get_three_prime_pair(int i){return mThreePrimePair.get_val(i);}
    const TS5Alignment<TRNAString>& get_alignment() {return mThreePrimePair.get_alignment();}

    // Method prototypes
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs,
            TS5SeedSites<TRNAString> const &pSeedSites);
    void clear_features();

private:
    TS5FeatSeedType<TRNAString> mSeedTypes;
    TS5FeatSitePos<TRNAString> mSitePos;
    TS5FeatAURich<TRNAString> mAURich;
    TS5FeatThreePrimePair<TRNAString> mThreePrimePair;
};

} // namespace ts5cs

#endif /* TS5_FEATURE_HPP_ */

