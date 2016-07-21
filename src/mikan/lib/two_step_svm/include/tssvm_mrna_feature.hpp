#ifndef TSSVM_MRNA_FEATURE_HPP_
#define TSSVM_MRNA_FEATURE_HPP_

#include <mikan/lib/two_step_svm/include/tssvm_seed_site.hpp>      // TSSVMSeedSites
#include <mikan/lib/two_step_svm/include/tssvm_site_svm.hpp>       // TSSVMSiteInputVector
#include <seqan/sequence.h>

namespace tssvm{

//
// UTR length feature
//
template <class TRNAString>
class TSSVMFeatUTRLen
{
public:
    // Define methods
    TSSVMFeatUTRLen() {}
    seqan::String<float>& get_val(int i){return mUTRLen[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mUTRLen;}
    void resize_feat(unsigned pLen){seqan::resize(mUTRLen, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mUTRLen;
};

//
// Site number feature
//
template <class TRNAString>
class TSSVMFeatSiteNum
{
public:
    // Define methods
    TSSVMFeatSiteNum() {}
    seqan::String<float>& get_val(int i){return mSiteNum[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSiteNum;}
    void resize_feat(unsigned pLen){seqan::resize(mSiteNum, pLen);seqan::resize(mSiteNumRaw, pLen);}
    seqan::String<unsigned>& get_site_count(){return mSiteNumRaw;}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mSiteNum;
    seqan::String<unsigned> mSiteNumRaw;
};

//
// Total discriminant / UTR length feature
//
template <class TRNAString>
class TSSVMFeatTotDiscUTRLen
{
public:
    // Define methods
    TSSVMFeatTotDiscUTRLen() {}
    seqan::String<float>& get_val(int i){return mTotDiscUTRLen[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mTotDiscUTRLen;}
    void resize_feat(unsigned pLen){seqan::resize(mTotDiscUTRLen, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mTotDiscUTRLen;

};

//
// Seed type number feature
//
template <class TRNAString>
class TSSVMFeatSeedTypeNum
{
public:
    // Define methods
    TSSVMFeatSeedTypeNum() {}
    seqan::String<float>& get_val(int i){return mSeedTypeNum[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSeedTypeNum;}
    void resize_feat(unsigned pLen){seqan::resize(mSeedTypeNum, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    seqan::StringSet<seqan::String<float> > mSeedTypeNum;

};

//
// Discriminant bins feature
//
template <class TRNAString>
class TSSVMFeatDiscBin
{
public:
    // Define methods
    TSSVMFeatDiscBin() {}
    seqan::String<float>& get_val(int i){return mDiscBin[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mDiscBin;}
    void resize_feat(unsigned pLen){seqan::resize(mDiscBin, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mDiscBin;

};

//
// Optimal distance feature
//
template <class TRNAString>
class TSSVMFeatOptDist
{
public:
    // Define methods
    TSSVMFeatOptDist() {}
    seqan::String<float>& get_val(int i){return mOptDist[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mOptDist;}
    void resize_feat(unsigned pLen){seqan::resize(mOptDist, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mOptDist;

};

//
// Site number flag feature
//
template <class TRNAString>
class TSSVMFeatSiteNumFlg
{
public:
    // Define methods
    TSSVMFeatSiteNumFlg() {}
    seqan::String<float>& get_val(int i){return mSiteNumFlg[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mSiteNumFlg;}
    void resize_feat(unsigned pLen){seqan::resize(mSiteNumFlg, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mSiteNumFlg;

};

//
// Total discriminant feature
//
template <class TRNAString>
class TSSVMFeatTotDisc
{
public:
    // Define methods
    TSSVMFeatTotDisc() {}
    seqan::String<float>& get_val(int i){return mTotDisc[i];}
    seqan::StringSet<seqan::String<float> >& get_all_val(){return mTotDisc;}
    void resize_feat(unsigned pLen){seqan::resize(mTotDisc, pLen);}

    // Method prototype
    int add_features(unsigned pMRNAPosIdx, seqan::String<unsigned>& pSitePosByMRNA,
            TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();

private:
    // Define variables
    seqan::StringSet<seqan::String<float> > mTotDisc;

};


//
// Store all raw feature values
//
template <class TRNAString>
class TSSVMRNARawFeatures
{
public:
    // Define types
    typedef seqan::StringSet<seqan::String<float> > TFeatSet;

    // Define variables
    seqan::String<bool> mEffectiveRNAs;

public:
    // Define methods
    TSSVMRNARawFeatures() {}
    seqan::String<float>& get_utr_len(int i){return mUTRLen.get_val(i);}
    seqan::String<float>& get_site_num(int i){return mSiteNum.get_val(i);}
    seqan::String<float>& get_tot_disc_utr_len(int i){return mTotDiscUTRLen.get_val(i);}
    seqan::String<float>& get_seed_type_num(int i){return mSeedTypeNum.get_val(i);}
    seqan::String<float>& get_disc_num(int i){return mDiscBin.get_val(i);}
    seqan::String<float>& get_opt_dist(int i){return mOptDist.get_val(i);}
    seqan::String<float>& get_site_num_flag(int i){return mSiteNumFlg.get_val(i);}
    seqan::String<float>& get_to_disc(int i){return mTotDisc.get_val(i);}

    TFeatSet& get_all_utr_len(){return mUTRLen.get_all_val();}
    TFeatSet& get_all_site_num(){return mSiteNum.get_all_val();}
    TFeatSet& get_all_tot_disc_utr_len(){return mTotDiscUTRLen.get_all_val();}
    TFeatSet& get_all_seed_type_num(){return mSeedTypeNum.get_all_val();}
    TFeatSet& get_all_disc_num(){return mDiscBin.get_all_val();}
    TFeatSet& get_all_opt_dist(){return mOptDist.get_all_val();}
    TFeatSet& get_all_site_num_flag(){return mSiteNumFlg.get_all_val();}
    TFeatSet& get_all_to_disc(){return mTotDisc.get_all_val();}

    seqan::String<unsigned>& get_site_count(){return mSiteNum.get_site_count();}

    // Method prototypes
    int add_features(TSSVMSeedSites<TRNAString> &pSeedSites, seqan::StringSet<TRNAString> const &pMRNASeqs,
            TSSVMSeedSiteOverlap<TRNAString> &pOverlappedSites,
            TSSVMSiteInputVector<TRNAString> &pSiteInput);
    void clear_features();
    void resize_feat(unsigned pLen);

private:
    typedef std::set<unsigned>::iterator TItSet;

    TSSVMFeatUTRLen<TRNAString> mUTRLen;
    TSSVMFeatSiteNum<TRNAString> mSiteNum;
    TSSVMFeatTotDiscUTRLen<TRNAString> mTotDiscUTRLen;
    TSSVMFeatSeedTypeNum<TRNAString> mSeedTypeNum;
    TSSVMFeatDiscBin<TRNAString> mDiscBin;
    TSSVMFeatOptDist<TRNAString> mOptDist;
    TSSVMFeatSiteNumFlg<TRNAString> mSiteNumFlg;
    TSSVMFeatTotDisc<TRNAString> mTotDisc;
};

} // namespace tssvm

#endif /* TSSVM_MRNA_FEATURE_HPP_ */
