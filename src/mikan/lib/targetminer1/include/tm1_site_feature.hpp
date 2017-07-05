#ifndef TM1_SITE_FEATURE_HPP_
#define TM1_SITE_FEATURE_HPP_

#include <tm1_align.hpp>         // TM1Alignment
#include <tm1_seed_site.hpp>     // TM1SeedSites
#include <tm1_site_cluster.hpp>  // TM1SortedSitePos
#include <seqan/sequence.h>

namespace tm1p {

//
// Seed type feature
//
template<class TRNAString>
class TM1FeatSeedType {
public:
    // Define types
    typedef typename seqan::StringSet<TRNAString> TRNASet;
    typedef typename seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TM1FeatSeedType() {}

    seqan::CharString &get_seed_type(unsigned idx) { return mSeedTypes[idx]; }

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites,
                     TM1SeedSites<TRNAString> &pSeedSites);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    seqan::StringSet<seqan::CharString> mSeedTypes;

private:
    void resize_features(unsigned pSize);
};

//
// Seed site position feature
//
template<class TRNAString>
class TM1FeatSitePos {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 15;
    static const int SEARCH_SEQ_LEN = 6;

    // Member variables
    seqan::String<int> mM1Pos;
    seqan::String<int> mM8Pos;
    seqan::String<int> mSeqLen;

public:
    // Define methods
    TM1FeatSitePos() {}

    int get_m1_pos(int i) { return mM1Pos[i]; }

    int get_m8_pos(int i) { return mM8Pos[i]; }

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites, TSitePos const &pSeedPos);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void resize_features(unsigned pSize);

};

//
// Distance to the nearest neighbor feature
//
template<class TRNAString>
class TM1FeatDistance {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TM1FeatDistance() {}

    int get_upstream_len(int i) { return mUpstream[i]; }

    int get_downstream_len(int i) { return mDownstream[i]; }

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1SortedSitePos<TRNAString> &pSortedSites);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    seqan::String<int> mUpstream;
    seqan::String<int> mDownstream;

private:
    void resize_features(unsigned pSize);

};

//
// AU-rich feature
//
template<class TRNAString>
class TM1FeatAURich {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<float> mAURich;

public:
    // Define methods
    TM1FeatAURich() {}

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1FeatDistance<TRNAString> &pDistance);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void calc_pos_scores(const TRNAString &pAURichRNA, int pStart, int pEnd, float &pTotalScore, float &pMaxScore);

    void resize_features(unsigned pSize);
};

//
// Single nucleotide frequency feature
//
template<class TRNAString>
class TM1FeatSingleFreq {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<int> mSeedCountU;
    seqan::String<int> mSeedCountC;
    seqan::String<int> mSeedTotal;

public:
    // Define methods
    TM1FeatSingleFreq() {}

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void resize_features(unsigned pSize);
};


//
// Single nucleotide frequency in flanking region feature
//
template<class TRNAString>
class TM1FeatSingleFreqFlank {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<int> m3PCountA;
    seqan::String<int> m3PCountU;
    seqan::String<int> m3PTotal;

public:
    // Define methods
    TM1FeatSingleFreqFlank() {}

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1FeatDistance<TRNAString> &pDistance);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void calc_pos_scores(const TRNAString &pMRNASeq, int pSiteStart, int pSiteEnd, int pIdx);

    void resize_features(unsigned pSize);
};

//
// Di-nucleotide frequency feature
//
template<class TRNAString>
class TM1FeatDiFreq {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<int> mSeedCountUC;
    seqan::String<int> mSeedCountCA;
    seqan::String<int> mSeedCountCG;
    seqan::String<int> mSeedTotal;

public:
    // Define methods
    TM1FeatDiFreq() {}

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void resize_features(unsigned pSize);
};

//
// Di-nucleotide frequency in flanking region feature
//
template<class TRNAString>
class TM1FeatDiFreqFlank {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<int> m3PCountAA;
    seqan::String<int> m3PCountAU;
    seqan::String<int> m3PCountAG;
    seqan::String<int> m3PCountAC;
    seqan::String<int> m3PCountUA;
    seqan::String<int> m3PCountUU;
    seqan::String<int> m3PCountUG;
    seqan::String<int> m3PCountUC;
    seqan::String<int> m3PCountGU;
    seqan::String<int> m3PCountCA;
    seqan::String<int> m3PCountCU;
    seqan::String<int> m3PTotal;

public:
    // Define methods
    TM1FeatDiFreqFlank() {}

    // Method prototype
    int add_features(TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites, TSitePos const &pMRNAPos,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1FeatDistance<TRNAString> &pDistance);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void calc_pos_scores(const TRNAString &pMRNASeq, int pSiteStart, int pSiteEnd, int pIdx);

    void resize_features(unsigned pSize);
};

//
// Single match frequency feature
//
template<class TRNAString>
class TM1FeatSingleMatch {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<int> mMatchUG;
    seqan::String<int> mMatchGU;
    seqan::String<int> mMatchCG;
    seqan::String<int> mMatchTotal;

public:
    // Define methods
    TM1FeatSingleMatch() {}

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites,
                     TSitePos const &pMRNAPos, TM1SeedSites<TRNAString> &pSeedSites,
                     TM1FeatSitePos<TRNAString> &pSeedPos);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void resize_features(unsigned pSize);

};

//
// Two consecutive matches frequency feature
//
template<class TRNAString>
class TM1FeatTwoConsecMatch {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

    // Member variables
    seqan::String<int> mMatchUACG;
    seqan::String<int> mMatchUAUG;
    seqan::String<int> mMatchCGGC;
    seqan::String<int> mMatchUGGC;
    seqan::String<int> mMatchTotal;

public:
    // Define methods
    TM1FeatTwoConsecMatch() {}

    // Method prototype
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs, seqan::String<bool> &pEffectiveSites,
                     TSitePos const &pMRNAPos, TM1SeedSites<TRNAString> &pSeedSites,
                     TM1FeatSitePos<TRNAString> &pSeedPos);

    void clear_features();

    void print_feature(unsigned pIdx);

private:
    void resize_features(unsigned pSize);
};

//
// Store all raw feature values
//
template<class TRNAString>
class TM1RawFeatures {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

public:
    // Define methods
    TM1RawFeatures() {}

    seqan::CharString &get_seed_type(int i) { return mSeedTypes.get_seed_type(i); }

    const TM1FeatAURich<TRNAString> &get_au_rich() { return mAURich; }

    const TM1FeatSingleFreq<TRNAString> &get_single_freq() { return mSingleFreqs; }

    const TM1FeatSingleFreqFlank<TRNAString> &get_single_freq_flank() { return mSingleFreqFlanks; }

    const TM1FeatDiFreq<TRNAString> &get_di_freq() { return mDiFreqs; }

    const TM1FeatDiFreqFlank<TRNAString> &get_di_freq_flank() { return mDiFreqFlanks; }

    const TM1FeatSingleMatch<TRNAString> &get_single_match() { return mSingleMatches; }

    const TM1FeatTwoConsecMatch<TRNAString> &get_two_consec_match() { return mTwoMatches; }

    const TM1Alignment<TRNAString> &get_alignment() { return mAlign; }

    // Method prototypes
    int add_features(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs,
                     TM1SeedSites<TRNAString> &pSeedSites, TM1SortedSitePos<TRNAString> &pSortedSites);

    void clear_features();

    void print_features(TM1SeedSites<TRNAString> &pSeedSites);

private:
    TM1FeatSeedType<TRNAString> mSeedTypes;
    TM1FeatSitePos<TRNAString> mSitePos;
    TM1FeatDistance<TRNAString> mDistance;
    TM1FeatAURich<TRNAString> mAURich;
    TM1FeatSingleFreq<TRNAString> mSingleFreqs;
    TM1FeatSingleFreqFlank<TRNAString> mSingleFreqFlanks;
    TM1FeatDiFreq<TRNAString> mDiFreqs;
    TM1FeatDiFreqFlank<TRNAString> mDiFreqFlanks;
    TM1FeatSingleMatch<TRNAString> mSingleMatches;
    TM1FeatTwoConsecMatch<TRNAString> mTwoMatches;

    TM1Alignment<TRNAString> mAlign;
};

} // namespace tm1p

#endif /* TM1_SITE_FEATURE_HPP_ */

