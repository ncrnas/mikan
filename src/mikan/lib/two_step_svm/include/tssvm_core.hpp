#ifndef TSSVM_CORE_HPP_
#define TSSVM_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"          // MKSequences
#include "mk_option.hpp"            // MKOptions
#include "tssvm_option.hpp"         // TSSVMOptions
#include "tssvm_seed_site.hpp"      // TSSVMSeedSites
#include "tssvm_site_filter.hpp"    // TSSVMSiteFilter
#include "tssvm_site_score.hpp"     // TSSVMSiteScores
#include "tssvm_rna_score.hpp"      // TSSVMRNAScores

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
    bool mExecSiteScore;
    bool mExecRNAFeat;
    bool mExecRNAScore;
    bool mOutputSiteScore;
    bool mOutputRNAScore;
    bool mOutputAlign;

    seqan::CharString mOFileSite;
    seqan::CharString mOFileRNA;

    mikan::TCharSet mSeedTypeDef;

public:
    // Define methods
    TSSVMCore(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
              mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
              mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecFilterOverlap(true),
            mExecSiteScore(true), mExecRNAFeat(true), mExecRNAScore(true), mOutputSiteScore(true),
            mOutputRNAScore(true), mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs),
            mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
            mRNAWithSites(pOpts), mSiteScores(pOpts), mSiteFilter(pOpts), mRNAScores(pOpts) {

        init_from_args(pOpts);

    }

    // Method prototypes
    void init_from_args(mikan::MKOptions const &opts);

    int open_output_file();

    int calculate_all_scores();

    int calculate_mirna_scores(unsigned pIdx);

private:
    mikan::TCharSet const &mMiRNAIds;
    mikan::TRNASet const &mMiRNASeqs;
    mikan::TCharSet const &mMRNAIds;
    mikan::TRNASet const &mMRNASeqs;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    TSSVMSeedSeqs mSeedSeqs;
    TSSVMSeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    TSSVMSiteScores mSiteScores;
    TSSVMSiteFilter mSiteFilter;
    TSSVMRNAScores mRNAScores;

private:
    int write_site_score(seqan::CharString const &pMiRNAId);

    int write_rna_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tssvm

#endif /* TSSVM_CORE_HPP_ */
