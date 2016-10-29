#ifndef TSSVM_SITE_FEATURE_HPP_
#define TSSVM_SITE_FEATURE_HPP_

#include <tssvm_align.hpp>          // TSAlign
#include <tssvm_seed_site.hpp>      // TSSVMSeedSite
#include <seqan/sequence.h>

namespace tssvm{

//
// Seed type feature
//
template <class TRNAString>
class TSSVMFeatSeedType
{
public:
    // Define methods
    TSSVMFeatSeedType() {}
    seqan::String<float>& get_val(int i){return mSeedTypes[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSeedTypes;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mSeedTypes;
};

//
// Alignment similarity feature
//
template <class TRNAString>
class TSSVMFeatSimilarity
{
public:
    // Define methods
    TSSVMFeatSimilarity() {}
    seqan::String<float>& get_val(int i){return mSimilarities[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSimilarities;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mSimilarities;
};

//
// 30nt upstream AU-rich feature
//
template <class TRNAString>
class TSSVMFeatAURichUp
{
public:
    // Define methods
    TSSVMFeatAURichUp() {}
    seqan::String<float>& get_val(int i){return mAURichUp[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mAURichUp;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mAURichUp;
    void getUpStreamPos(unsigned pStartMRNAPos, unsigned &pStartPos, unsigned &pEndPos);

};

//
// 30nt downstream AU-rich feature
//
template <class TRNAString>
class TSSVMFeatAURichDown
{
public:
    // Define methods
    TSSVMFeatAURichDown() {}
    seqan::String<float>& get_val(int i){return mAURichDown[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mAURichDown;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mAURichDown;
    void getDownStreamPos(unsigned pStartMRNAPos, unsigned pSeqLen, unsigned &pStartPos, unsigned &pEndPos);

};

//
// Seed site position feature
//
template <class TRNAString>
class TSSVMFeatSitePos
{
public:
    // Define methods
    TSSVMFeatSitePos() {}
    seqan::String<float>& get_val(int i){return mSitePos[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSitePos;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mSitePos;

};

//
// Sequence match feature
//
template <class TRNAString>
class TSSVMFeatSeqMatch
{
public:
    // Define methods
    TSSVMFeatSeqMatch() {}
    seqan::String<float>& get_val(int i){return mSeqMatch[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSeqMatch;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mSeqMatch;

};

//
// A1 match feature
//
template <class TRNAString>
class TSSVMFeatA1Match
{
public:
    // Define methods
    TSSVMFeatA1Match() {}
    seqan::String<float>& get_val(int i){return mA1Match[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mA1Match;}

    // Method prototype
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs,
            seqan::String<bool> &pEffectiveSites);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mA1Match;

};

//
// Store all raw feature values
//
template <class TRNAString>
class TSSVMRawFeatures
{
public:
    // Define types
    typedef seqan::StringSet<seqan::String<float> > TFeatSet;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    TSSVMRawFeatures() {}
    seqan::String<float>& get_seed_type(int i){return mSeedTypes.get_val(i);}
    seqan::String<float>& get_similarities(int i){return mSimilarities.get_val(i);}
    seqan::String<float>& get_au_rich_up(int i){return mAURichUp.get_val(i);}
    seqan::String<float>& get_au_rich_down(int i){return mAURichDown.get_val(i);}
    seqan::String<float>& get_site_pos(int i){return mSitePos.get_val(i);}
    seqan::String<float>& get_seq_match(int i){return mSeqMatch.get_val(i);}
    seqan::String<float>& get_a1_match(int i){return mA1Match.get_val(i);}

    TFeatSet& get_all_seed_type(){return mSeedTypes.get_all_val();}
    TFeatSet& get_all_similarities(){return mSimilarities.get_all_val();}
    TFeatSet& get_all_au_rich_up(){return mAURichUp.get_all_val();}
    TFeatSet& get_all_au_rich_down(){return mAURichDown.get_all_val();}
    TFeatSet& get_all_site_pos(){return mSitePos.get_all_val();}
    TFeatSet& get_all_seq_match(){return mSeqMatch.get_all_val();}
    TFeatSet& get_all_a1_match(){return mA1Match.get_all_val();}

    // Method prototypes
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, TSAlign<TRNAString> const &pAlignSeqs,
            TRNAString const &pMiRNASeq, seqan::StringSet<TRNAString> const &pMRNASeqs);
    void clear_features();

private:
    TSSVMFeatSeedType<TRNAString> mSeedTypes;
    TSSVMFeatSimilarity<TRNAString> mSimilarities;
    TSSVMFeatAURichUp<TRNAString> mAURichUp;
    TSSVMFeatAURichDown<TRNAString> mAURichDown;
    TSSVMFeatSitePos<TRNAString> mSitePos;
    TSSVMFeatSeqMatch<TRNAString> mSeqMatch;
    TSSVMFeatA1Match<TRNAString> mA1Match;
};

} // namespace tssvm

#endif /* TSSVM_SITE_FEATURE_HPP_ */
