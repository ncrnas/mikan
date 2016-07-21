#ifndef TS5_SCORE_HPP_
#define TS5_SCORE_HPP_

#include <mikan/lib/targetscan5_cs/include/ts5_feature.hpp>       // TS5RawFeatures
#include <mikan/lib/targetscan5_cs/include/ts5_seed_site.hpp>     // TS5SeedSites
#include <seqan/sequence.h>

namespace ts5cs{

//
// Seed type score
//
class TS5ScoreSeedType
{
public:
    // Define methods
    TS5ScoreSeedType(): mContrib8mer(-0.310), mContrib7M8(-0.161), mContrib7A1(-0.099) {}
    void set_score(int i, float val){mSeedTypes[i] = val;};
    float& get_score(int i){return mSeedTypes[i];}

    // Method prototypes
    void clear_scores();
    float calc_score(int pIdx, seqan::CharString &pSeedType);
    void resize_scores(int pSize);

private:
    seqan::String<float> mSeedTypes;
    float mContrib8mer;
    float mContrib7M8;
    float mContrib7A1;

};

//
// Seed site position score
//
class TS5ScoreSitePos
{
public:
    // Define methods
    TS5ScoreSitePos(): mSlope8mer(0.000172), mSlope7M8(0.000091), mSlope7A1(0.000072),
    mICept8mer(-0.07), mICept7M8(-0.037), mICept7A1(-0.032) {}
    void set_score(int i, float val){mSitePos[i] = val;};
    float& get_score(int i){return mSitePos[i];}

    // Method prototypes
    void clear_scores();
    float calc_score(int pIdx, seqan::CharString &pSeedType, int &pSitePos);
    void resize_scores(int pSize);

private:
    seqan::String<float> mSitePos;
    float mSlope8mer;
    float mSlope7M8;
    float mSlope7A1;
    float mICept8mer;
    float mICept7M8;
    float mICept7A1;
};

//
// AU-rich score
//
class TS5ScoreAURich
{
public:
    // Define methods
    TS5ScoreAURich(): mSlope8mer(-0.64), mSlope7M8(-0.5), mSlope7A1(-0.42),
    mICept8mer(0.365), mICept7M8(0.269), mICept7A1(0.236) {}
    void set_score(int i, float val){mAURich[i] = val;};
    float& get_score(int i){return mAURich[i];}

    // Method prototypes
    void clear_scores();
    float calc_score(int pIdx, seqan::CharString &pSeedType, float &pAuRich);
    void resize_scores(int pSize);

private:
    seqan::String<float> mAURich;
    float mSlope8mer;
    float mSlope7M8;
    float mSlope7A1;
    float mICept8mer;
    float mICept7M8;
    float mICept7A1;
};

//
//  Additional 3' paring score
//
class TS5ScoreThreePrimePair
{
public:
    // Define methods
    TS5ScoreThreePrimePair(): mSlope8mer(-0.0041), mSlope7M8(-0.031), mSlope7A1(-0.0211),
    mICept8mer(0.011), mICept7M8(0.067), mICept7A1(0.046) {}
    void set_score(int i, float val){mThreePrimePair[i] = val;}
    float& get_score(int i){return mThreePrimePair[i];}

    // Method prototypes
    void clear_scores();
    float calc_score(int pIdx, seqan::CharString &pSeedType, float &pThreePrimePair);
    void resize_scores(int pSize);

private:
    seqan::String<float> mThreePrimePair;
    float mSlope8mer;
    float mSlope7M8;
    float mSlope7A1;
    float mICept8mer;
    float mICept7M8;
    float mICept7A1;

};

//
// Store context scores
//
template <class TRNAString>
class TS5ContextScores
{
public:
    // Define methods
    TS5ContextScores() {}
    void set_score(int i, float val){mContextScores[i] = val;};

    float const& get_score(int i) const {return mContextScores[i];}
    float& get_seed_type_score(int i){return mSeedTypes.get_score(i);}
    float& get_site_pos_score(int i){return mSitePos.get_score(i);}
    float& get_au_rich_score(int i){return mAURich.get_score(i);}
    float& get_three_prime_pair_score(int i){return mThreePrimePair.get_score(i);}

    // Define variables
    seqan::String<bool> mEffectiveSites;

    // Method prototypes
    void clear_scores();
    int calc_scores(TS5RawFeatures<TRNAString> &pRawFeatures);
    void resize_scores(int pSize);

private:
    seqan::String<float> mContextScores;
    TS5ScoreSeedType mSeedTypes;
    TS5ScoreSitePos mSitePos;
    TS5ScoreAURich mAURich;
    TS5ScoreThreePrimePair mThreePrimePair;

};

//
// Total context scores
//
template <class TRNAString>
class TS5TotalScores
{
public:
    // Define methods
    TS5TotalScores() {}
    const seqan::String<float>& get_scores(){return mTotalScores;}
    const seqan::String<int>& get_mrna_pos(){return mMRNAPos;}
    const seqan::String<int>& get_site_num(){return mSiteNum;}

    // Method prototypes
    void clear_scores();
    int calc_scores(TS5SeedSites<TRNAString> &pSeedSites, TS5ContextScores<TRNAString> &pContextScores);

private:
    seqan::String<float> mTotalScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace ts5cs

#endif /* TS5_SCORE_HPP_ */