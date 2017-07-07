#ifndef TSSVM_CORE_HPP_
#define TSSVM_CORE_HPP_

#include "tssvm_align.hpp"          // TSAlign
#include "mk_typedef.hpp"           // TRNATYPE
#include "tssvm_mrna_feature.hpp"   // TSSVMRNARawFeatures
#include "tssvm_mrna_svm.hpp"       // TSSVMRNAInputVector
#include "tssvm_option.hpp"         // TSSVMOptions
#include "tssvm_seed_site.hpp"      // TSSVMSeedSites, TSSVMSeedSiteOverlap
#include "tssvm_site_svm.hpp"       // TSSVMSiteInputVector
#include "mk_sequence.hpp"          // MKSequences
#include "mk_typedef.hpp"           // TCharSet, TRNASet, TIndexQGram, TFinder
#include <seqan/sequence.h>

namespace tssvm {

int TSSVMCoreMain(int argc, char const **argv);

//
//Two-step SVM score process core
//
class TSSVMCore {
public:
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

    seqan::CharString mOFileTargetSite;
    seqan::CharString mOFileMRNA;

public:
    // Define methods
    TSSVMCore(mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs, 
              mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
              mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
              mExecSearchSeedSites(true), mExecFilterOverlap(true), mExecAlignSeq(true), mExecSiteFeat(true),
              mExecSiteScore(true), mExecRNAFeat(true), mExecRNAScore(true), mOutputSiteScore(true),
              mOutputRNAScore(true), mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs),
              mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
              mSiteInput(mSiteModel) {}

    // Method prototypes
    void init_from_args(TSSVMOptions &opts);

    int open_output_file();

    int calculate_all_scores();

    int calculate_mirna_scores(unsigned pIdx);

    int init_site_svm();

private:
    mikan::TCharSet const &mMiRNAIds;
    mikan::TRNASet const &mMiRNASeqs;
    mikan::TCharSet const &mMRNAIds;
    mikan::TRNASet const &mMRNASeqs;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    TSSVMSeedSites<mikan::TRNATYPE> mSeedSites;
    TSSVMSeedSiteOverlap<mikan::TRNATYPE> mOverlappedSites;
    TSAlign<mikan::TRNATYPE> mAlignSeqs;
    TSSVMRawFeatures<mikan::TRNATYPE> mSiteFeatures;
    TSSVMSiteModel<mikan::TRNATYPE> mSiteModel;
    TSSVMSiteInputVector<mikan::TRNATYPE> mSiteInput;
    TSSVMRNARawFeatures<mikan::TRNATYPE> mRnaFeatures;
    TSSVMRNAInputVector<mikan::TRNATYPE> mRnaInput;

private:
    int write_ts_scores(seqan::CharString const &pMiRNAId);

    int write_mrna_scores(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tssvm

#endif /* TSSVM_CORE_HPP_ */
