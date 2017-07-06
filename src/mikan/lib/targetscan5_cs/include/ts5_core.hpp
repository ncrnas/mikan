#ifndef TS5_CORE_HPP_
#define TS5_CORE_HPP_

#include <ts5_feature.hpp>       // TS5RawFeatures
#include <ts5_inst_template.hpp> // TRNATYPE
#include <ts5_option.hpp>        // TS5CSOptions
#include <ts5_score.hpp>         // TS5ContextScores, TS5TotalScores
#include <ts5_seed_site.hpp>     // TS5SeedSites
#include <mk_sequence.hpp>       // MKSequences

namespace ts5cs {

int TS5CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
template<class TRNAString, int SEEDLEN = 6>
class TS5Core {
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecGetRawFeat;
    bool mExecCalcContexScore;
    bool mExecSumScores;
    bool mOutputContexScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileContext;
    seqan::CharString mOFileTotal;

public:
    // Define methods
    TS5Core(TCharSet const &pMiRNAIds, TRNASet const &pMiRNASeqs, TCharSet const &pMRNAIds,
            TRNASet const &pMRNASeqs,
            TIndexQGram &pRNAIdx, TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecGetRawFeat(true), mExecCalcContexScore(true),
            mExecSumScores(true), mOutputContexScore(true), mOutputTotalScore(true),
            mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
            mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder) {}

    // Method prototypes
    void init_from_args(TS5CSOptions &opts);

    int open_output_file();

    int calculate_all_scores();

    int calculate_mirna_scores(unsigned pIdx);

private:
    TCharSet const &mMiRNAIds;
    TRNASet const &mMiRNASeqs;
    TCharSet const &mMRNAIds;
    TRNASet const &mMRNASeqs;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    TS5SeedSites<TRNAString> mSeedSites;
    TS5RawFeatures<TRNAString> mRawFeatures;
    TS5ContextScores<TRNAString> mCsScores;
    TS5TotalScores<TRNAString> mTotalScore;

private:
    int write_context_score(seqan::CharString const &pMiRNAId);

    int write_total_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace ts5cs

#endif /* TS5_CORE_HPP_ */
