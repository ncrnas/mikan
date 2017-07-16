#ifndef TM1_CORE_HPP_
#define TM1_CORE_HPP_

#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "tm1_mrna_feature.hpp"  // TM1MRNAFeatures
#include "tm1_mrna_svm.hpp"      // TM1MRNAModel, TM1MRNAInputVector
#include "tm1_option.hpp"        // TM1CSOptions
#include "tm1_score.hpp"         // TM1ClassifiedScores
#include "tm1_seed_site.hpp"     // TM1SeedSites
#include "tm1_site_cluster.hpp"  // TM1Overlap
#include "tm1_site_feature.hpp"  // TM1RawFeatures

namespace tm1p {

int TM1CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
class TM1Core {
public:
    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecGetRawFeat;
    bool mExecSortSites;
    bool mExecGetMRNAFeat;
    bool mExecRNAScore;
    bool mExecSumScores;
    bool mOutputSitePos;
    bool mOutputScore;
    bool mOutputAlign;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileScore;

    mikan::TCharSet mSeedTypeDef;

public:
    // Define methods
    TM1Core(mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs, 
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecGetRawFeat(true), mExecSortSites(true), mExecGetMRNAFeat(true),
            mExecRNAScore(true), mExecSumScores(true), mOutputSitePos(true), mOutputScore(true), mOutputAlign(true),
            mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
            mSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

    // Method prototypes
    void init_from_args(TM1CSOptions &opts);

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

    TM1SeedSites mSeedSites;
    TM1RawFeatures mRawFeatures;
    TM1SortedSitePos mSortedSites;
    TM1MRNAFeatures mMRNAFeatures;
    TM1MRNAInputVector mMRNAInput;
    TM1ClassifiedScores mScores;

private:
    int write_site_positions(seqan::CharString const &pMiRNAId);

    int write_scores(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tm1p

#endif /* TM1_CORE_HPP_ */
