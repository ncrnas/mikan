#ifndef MR3_CORE_HPP_
#define MR3_CORE_HPP_

#include <mr3_inst_template.hpp>  // TRNATYPE
#include <mr3_option.hpp>         // MR3Options
#include <mr3_score.hpp>          // MR3GGDScores, MR3TotalScores
#include <mr3_seed_site.hpp>      // MR3Sequences, MR3SeedSites
#include <mr3_site_cluster.hpp>   // MR3Overlap, MR3SortedSitePos
#include <seqan/sequence.h>

namespace mr3as {

//
// Input data for MR3 score
//
template <class TRNAString>
class MR3CoreInput
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
    MR3CoreInput() {}
    TCharSet const& get_mirna_ids() {return mMiRNASeqs.get_ids();}
    TRNASet const& get_mirna_seqs () {return mMiRNASeqs.get_seqs();}
    TCharSet const& get_mrna_ids() {return mMRNASeqs.get_ids();}
    TRNASet const& get_mrna_seqs () {return mMRNASeqs.get_seqs();}

    // Method prototypes
    void init_from_args(MR3Options& opts);
    int load_seq_from_file();

private:
    MR3Sequences<TRNAString> mMiRNASeqs;
    MR3Sequences<TRNAString> mMRNASeqs;
};

//
// MR3 score process core
//
template <class TRNAString, int SEEDLEN=6>
class MR3Core
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Constant values
    static const unsigned INDEXED_SEQ_LEN = SEEDLEN;
    static const unsigned OVERLAP_LEN = 6;

    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalSiteScore;
    bool mExecFilterOverlap;
    bool mExecSortSites;
    bool mExecSumScores;
    bool mOutputSiteScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileTotal;
    int mMinSeedLen;
    int mMaxSeedLen;
    float mMinAlignScore;
    float mMaxEnergy;

    seqan::StringSet<seqan::CharString> mSeedDef;

public:
    // Define methods
    MR3Core(TCharSet const& pMiRNAIds, TRNASet const& pMiRNASeqs, TCharSet const& pMRNAIds, TRNASet const& pMRNASeqs,
            TIndexQGram& pRNAIdx, TFinder& pFinder):
                mExecSearchSeedSites(true), mExecCalSiteScore(true), mExecFilterOverlap(true),
                mExecSortSites(true), mExecSumScores(true), mOutputSiteScore(true), mOutputTotalScore(true),
                mOutputAlign(true), mMinSeedLen(6), mMaxSeedLen(8), mMinAlignScore(120.0), mMaxEnergy(1.0),
                mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
                mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
                mSiteScores()
                {}

    // Method prototypes
    void init_from_args(MR3Options& opts);
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

    MR3SeedSites<TRNAString> mSeedSites;
    MR3SiteScores<TRNAString> mSiteScores;
    MR3Overlap<TRNAString> mOverlappedSites;
    MR3SortedSitePos<TRNAString> mSortedSites;
    MR3TotalScores<TRNAString> mTotalScores;

private:
    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_total_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace mr3as

#endif /* MR3_CORE_HPP_ */
