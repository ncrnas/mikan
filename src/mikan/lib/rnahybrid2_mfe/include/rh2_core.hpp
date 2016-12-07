#ifndef RH2_CORE_HPP_
#define RH2_CORE_HPP_

#include <rh2_inst_template.hpp> // TRNATYPE
#include <rh2_option.hpp>        // RH2Options
#include <rh2_score.hpp>         // RH2MFEScores, RH2TotalScores
#include <rh2_seed_site.hpp>     // RH2Sequences, RH2SeedSites
#include <rh2_site_cluster.hpp>  // RH2Overlap, RH2TopNScore, RH2SortedSitePos
#include <seqan/sequence.h>

namespace rh2mfe {

int RH2CoreMain(int argc, char const ** argv);

//
// Input data for RNAhybrid MFE score
//
template <class TRNAString>
class RH2CoreInput
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
    RH2CoreInput() {}
    TCharSet const& get_mirna_ids() {return mMiRNASeqs.get_ids();}
    TRNASet const& get_mirna_seqs () {return mMiRNASeqs.get_seqs();}
    TCharSet const& get_mrna_ids() {return mMRNASeqs.get_ids();}
    TRNASet const& get_mrna_seqs () {return mMRNASeqs.get_seqs();}

    // Method prototypes
    void init_from_args(RH2Options& opts);
    int load_seq_from_file();

private:
    RH2Sequences<TRNAString> mMiRNASeqs;
    RH2Sequences<TRNAString> mMRNASeqs;
};

//
// RNAhybrid MFE score process core
//
template <class TRNAString, int SEEDLEN=6>
class RH2Core
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalMFEScore;
    bool mExecFilterOverlap;
    bool mExecFilterSiteNum;
    bool mExecSortSites;
    bool mExecSumScores;
    bool mOutputMFEScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileMFE;
    seqan::CharString mOFileTotal;

    seqan::CharString mSeedDef;
    seqan::CharString mOverlapDef;
    int mMaxHits;

public:
    // Define methods
    RH2Core(TCharSet const& pMiRNAIds, TRNASet const& pMiRNASeqs, TCharSet const& pMRNAIds, TRNASet const& pMRNASeqs,
            TIndexQGram& pRNAIdx, TFinder& pFinder, int pMRNAMaxLen, int pMiRNAMaxLen, std::string& pSeedDef):
                mExecSearchSeedSites(true), mExecCalMFEScore(true), mExecFilterOverlap(true),
                mExecFilterSiteNum(true), mExecSortSites(true), mExecSumScores(true), mOutputMFEScore(true),
                mOutputTotalScore(true), mOutputAlign(true), mMaxHits(0), mMiRNAIds(pMiRNAIds),
                mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
                mSeedSites(pRNAIdx, pFinder, pMRNASeqs), mMfeScores(pMRNAMaxLen, pMiRNAMaxLen, pSeedDef)
                {}

    // Method prototypes
    void init_from_args(RH2Options& opts);
    int open_output_file();
    int calculate_all_scores();
    int calculate_mirna_scores(unsigned pIdx);

private:
    TCharSet const& mMiRNAIds;
    TRNASet const& mMiRNASeqs;
    TCharSet const& mMRNAIds;
    TRNASet const& mMRNASeqs;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    RH2SeedSites<TRNAString> mSeedSites;
    RH2MFEScores<TRNAString> mMfeScores;
    RH2Overlap<TRNAString> mOverlappedSites;
    RH2TopNScore<TRNAString> mTopScoredSites;
    RH2SortedSitePos<TRNAString> mSortedSites;
    RH2TotalScores<TRNAString> mTotalScores;

private:
    int write_mfe_score(seqan::CharString const &pMiRNAId);
    int write_total_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace rh2mfe

#endif /* RH2_CORE_HPP_ */
