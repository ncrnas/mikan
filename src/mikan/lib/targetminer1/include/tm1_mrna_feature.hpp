#ifndef TM1_MRNA_FEATURE_HPP_
#define TM1_MRNA_FEATURE_HPP_

#include <tm1_seed_site.hpp>     // TM1SeedSites
#include <tm1_site_cluster.hpp>  // TM1SortedSitePos
#include <tm1_site_feature.hpp>  // TM1RawFeatures
#include <seqan/sequence.h>

namespace tm1p {

//
// Seed type feature
//
template<class TRNAString>
class TM1MRNASeedType {
public:
    // Member variables
    seqan::String<float> mNum6mer;
    seqan::String<float> mNum8mer;

public:
    // Define methods
    TM1MRNASeedType() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};

//
// Site count feature
//
template<class TRNAString>
class TM1MRNASiteCount {
public:
    // Member variables
    seqan::String<unsigned> mSiteCount;

public:
    // Define methods
    TM1MRNASiteCount() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};

//
// AU-rich feature
//
template<class TRNAString>
class TM1MRNAAURich {
public:
    seqan::String<float> mAU6mer;
    seqan::String<float> mAU7merm8;
    seqan::String<float> mAU8mer;

public:
    // Define methods
    TM1MRNAAURich() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};

//
// Single nucleotide frequency feature
//
template<class TRNAString>
class TM1MRNASingleFreq {
public:
    seqan::String<float> mSeedFreqU;
    seqan::String<float> mSeedFreqC;

public:
    // Define methods
    TM1MRNASingleFreq() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};


//
// Single nucleotide frequency in flanking region feature
//
template<class TRNAString>
class TM1MRNASingleFreqFlank {
public:
    seqan::String<float> m3PFreqA;
    seqan::String<float> m3PFreqU;

public:
    // Define methods
    TM1MRNASingleFreqFlank() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};

//
// Di-nucleotide frequency feature
//
template<class TRNAString>
class TM1MRNADiFreq {
public:
    seqan::String<float> mSeedFreqUC;
    seqan::String<float> mSeedFreqCA;
    seqan::String<float> mSeedFreqCG;

public:
    // Define methods
    TM1MRNADiFreq() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};

//
// Di-nucleotide frequency in flanking region feature
//
template<class TRNAString>
class TM1MRNADiFreqFlank {
public:
    seqan::String<float> m3PFreqAA;
    seqan::String<float> m3PFreqAU;
    seqan::String<float> m3PFreqAG;
    seqan::String<float> m3PFreqAC;
    seqan::String<float> m3PFreqUA;
    seqan::String<float> m3PFreqUU;
    seqan::String<float> m3PFreqUG;
    seqan::String<float> m3PFreqUC;
    seqan::String<float> m3PFreqGU;
    seqan::String<float> m3PFreqCA;
    seqan::String<float> m3PFreqCU;

public:
    // Define methods
    TM1MRNADiFreqFlank() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};

//
// Single match frequency feature
//
template<class TRNAString>
class TM1MRNASingleMatch {
public:
    seqan::String<float> mFreqUG;
    seqan::String<float> mFreqGU;
    seqan::String<float> mFreqCG;

public:
    // Define methods
    TM1MRNASingleMatch() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);
};

//
// Two consecutive matches frequency feature
//
template<class TRNAString>
class TM1MRNATwoConsecMatch {
public:
    seqan::String<float> mFreqUACG;
    seqan::String<float> mFreqUAUG;
    seqan::String<float> mFreqCGGC;
    seqan::String<float> mFreqUGGC;

public:
    // Define methods
    TM1MRNATwoConsecMatch() {}

    // Method prototype
    int add_features(unsigned pIdx, const seqan::String<unsigned> &pSortedSites,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_feature(unsigned pIdx);

};


//
// Scale and store all feature values
//
template<class TRNAString>
class TM1ScaledFeatures {
public:
    // Constant values
    static const unsigned FEATURE_NUM = 30;

public:
    // Define methods
    TM1ScaledFeatures() : mLower(-1.0), mUpper(1.0) {
        init_max_vals();
    }

    const seqan::String<unsigned> &get_mrna_ids() { return mMRNAIDs; }

    const seqan::StringSet<seqan::String<float> > &get_mrna_features() { return mScaledFeatures; }

    // Method prototypes
    int scale_features(TM1SortedSitePos<TRNAString> &pSortedSites, TM1MRNASeedType<TRNAString> &pSeedTypes,
                       TM1MRNAAURich<TRNAString> &pAURich, TM1MRNASingleFreq<TRNAString> &pSingleFreqs,
                       TM1MRNASingleFreqFlank<TRNAString> &pSingleFreqFlanks, TM1MRNADiFreq<TRNAString> &pDiFreqs,
                       TM1MRNADiFreqFlank<TRNAString> &pDiFreqFlanks,
                       TM1MRNASingleMatch<TRNAString> &pSingleMatches,
                       TM1MRNATwoConsecMatch<TRNAString> &pTwoMatches);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_features();

private:
    seqan::String<unsigned> mMRNAIDs;
    seqan::StringSet<seqan::String<float> > mScaledFeatures;

    seqan::String<float> mMaxFeatVals;
    const float mLower;
    const float mUpper;

private:
    void init_max_vals();

    float scale_feat_val(unsigned pIdx, float pVal);

    void scale_seed_type(unsigned pIdx, TM1MRNASeedType<TRNAString> &pSeedTypes);

    void scale_au_rich(unsigned pIdx, TM1MRNAAURich<TRNAString> &pAURich);

    void scale_single_freq(unsigned pIdx, TM1MRNASingleFreq<TRNAString> &pSingleFreqs);

    void scale_freq_flank(unsigned pIdx, TM1MRNASingleFreqFlank<TRNAString> &pSingleFreqFlanks);

    void scale_di_freq(unsigned pIdx, TM1MRNADiFreq<TRNAString> &pDiFreqs);

    void scale_di_freq_flank(unsigned pIdx, TM1MRNADiFreqFlank<TRNAString> &pDiFreqFlanks);

    void scale_single_match(unsigned pIdx, TM1MRNASingleMatch<TRNAString> &pSingleMatches);

    void scale_two_match(unsigned pIdx, TM1MRNATwoConsecMatch<TRNAString> &pTwoMatches);
};

//
// Store all raw feature values
//
template<class TRNAString>
class TM1MRNAFeatures {
public:
    // Define methods
    TM1MRNAFeatures() {}

    TM1ScaledFeatures<TRNAString> &get_scaled_feature() { return mScaledFeats; }

    const seqan::String<unsigned> &get_site_counts() { return mSiteCounts.mSiteCount; }

    // Method prototypes
    int add_features(TM1SeedSites<TRNAString> &pSeedSites, TM1RawFeatures<TRNAString> &pRawFeatures,
                     TM1SortedSitePos<TRNAString> &pSortedSites);

    void clear_features();

    void resize_features(unsigned pSize);

    void print_features(TM1SortedSitePos<TRNAString> &pSortedSites);

private:
    TM1MRNASeedType<TRNAString> mSeedTypes;
    TM1MRNASiteCount<TRNAString> mSiteCounts;
    TM1MRNAAURich<TRNAString> mAURich;
    TM1MRNASingleFreq<TRNAString> mSingleFreqs;
    TM1MRNASingleFreqFlank<TRNAString> mSingleFreqFlanks;
    TM1MRNADiFreq<TRNAString> mDiFreqs;
    TM1MRNADiFreqFlank<TRNAString> mDiFreqFlanks;
    TM1MRNASingleMatch<TRNAString> mSingleMatches;
    TM1MRNATwoConsecMatch<TRNAString> mTwoMatches;

    TM1ScaledFeatures<TRNAString> mScaledFeats;
};

} // namespace tm1p

#endif /* TM1_MRNA_FEATURE_HPP_ */
