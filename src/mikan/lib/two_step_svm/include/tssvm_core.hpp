#ifndef TSSVM_CORE_HPP_
#define TSSVM_CORE_HPP_

#include <mikan/lib/two_step_svm/include/tssvm_align.hpp>          // TSAlign
#include <mikan/lib/two_step_svm/include/tssvm_mrna_feature.hpp>   // TSSVMRNARawFeatures
#include <mikan/lib/two_step_svm/include/tssvm_mrna_svm.hpp>       // TSSVMRNAInputVector
#include <mikan/lib/two_step_svm/include/tssvm_option.hpp>         // TSSVMOptions
#include <mikan/lib/two_step_svm/include/tssvm_seed_site.hpp>      // TSSVMSequences, TSSVMSeedSites, TSSVMSeedSiteOverlap
#include <mikan/lib/two_step_svm/include/tssvm_site_svm.hpp>       // TSSVMSiteInputVector
#include <seqan/sequence.h>

namespace tssvm{

//
// Input data for Two-step SVM score
//
template <class TRNAString>
class TSSVMCoreInput
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Declare variables
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;

public:
    // Define methods
    TSSVMCoreInput() {}
    TCharSet const& get_mirna_ids() {return mMiRNASeqs.get_ids();}
    TRNASet const& get_mirna_seqs () {return mMiRNASeqs.get_seqs();}
    TCharSet const& get_mrna_ids() {return mMRNASeqs.get_ids();}
    TRNASet const& get_mrna_seqs () {return mMRNASeqs.get_seqs();}

    // Method prototypes
    void init_from_args(TSSVMOptions& opts);
    int load_seq_from_file();

private:
    TSSVMSequences<TRNAString> mMiRNASeqs;
    TSSVMSequences<TRNAString> mMRNASeqs;
};

//
//Two-step SVM score process core
//
template <class TRNAString, int SEEDLEN=6>
class TSSVMCore
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecFilterOverlap;
    bool mExecAlignSeq;
    bool mExecSiteFeat;
    bool mExecSiteScore;
    bool mExecRNAFeat;
    bool mExecRNAScore;
    bool mOutputSiteScore;
    bool mOutputRNAScore;
    bool mOutputAlign;

    seqan::CharString mModelPath;
    seqan::CharString mOFileTargetSite;
    seqan::CharString mOFileMRNA;

public:
    // Define methods
    TSSVMCore(TCharSet const& pMiRNAIds, TRNASet const& pMiRNASeqs, TCharSet const& pMRNAIds,
            TRNASet const& pMRNASeqs, TIndexQGram& pRNAIdx, TFinder& pFinder):
                  mExecSearchSeedSites(true), mExecFilterOverlap(true), mExecAlignSeq(true), mExecSiteFeat(true),
                  mExecSiteScore(true), mExecRNAFeat(true), mExecRNAScore(true), mOutputSiteScore(true),
                  mOutputRNAScore(true), mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs),
                  mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
                  mSiteInput(mSiteModel) {}

    // Method prototypes
    void init_from_args(TSSVMOptions& opts);
    int open_output_file();
    int calculate_all_scores();
    int calculate_mirna_scores(unsigned pIdx);
    int init_site_svm();

private:
    TCharSet const& mMiRNAIds;
    TRNASet const& mMiRNASeqs;
    TCharSet const& mMRNAIds;
    TRNASet const& mMRNASeqs;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    TSSVMSeedSites<TRNAString> mSeedSites;
    TSSVMSeedSiteOverlap<TRNAString> mOverlappedSites;
    TSAlign<TRNAString> mAlignSeqs;
    TSSVMRawFeatures<TRNAString> mSiteFeatures;
    TSSVMSiteModel<TRNAString> mSiteModel;
    TSSVMSiteInputVector<TRNAString> mSiteInput;
    TSSVMRNARawFeatures<TRNAString> mRnaFeatures;
    TSSVMRNAInputVector<TRNAString> mRnaInput;

private:
    int write_ts_scores(seqan::CharString const &pMiRNAId);
    int write_mrna_scores(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tssvm

#endif /* TSSVM_CORE_HPP_ */