#ifndef TM1_CORE_HPP_
#define TM1_CORE_HPP_

#include <tm1_inst_template.hpp> // TRNATYPE
#include <tm1_mrna_feature.hpp>  // TM1MRNAFeatures
#include <tm1_mrna_svm.hpp>      // TM1MRNAModel, TM1MRNAInputVector
#include <tm1_option.hpp>        // TM1CSOptions
#include <tm1_score.hpp>         // TM1ClassifiedScores
#include <tm1_seed_site.hpp>     // TM1Sequences, TM1SeedSites
#include <tm1_site_cluster.hpp>  // TM1Overlap
#include <tm1_site_feature.hpp>  // TM1RawFeatures

namespace tm1p{

int TM1CoreMain(int argc, char const ** argv);

//
// Input data for TargetScan context score
//
template <class TRNAString>
class TM1CoreInput
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
    TM1CoreInput() {}
    TCharSet const& get_mirna_ids() {return mMiRNASeqs.get_ids();}
    TRNASet const& get_mirna_seqs () {return mMiRNASeqs.get_seqs();}
    TCharSet const& get_mrna_ids() {return mMRNASeqs.get_ids();}
    TRNASet const& get_mrna_seqs () {return mMRNASeqs.get_seqs();}

    // Method prototypes
    void init_from_args(TM1CSOptions& opts);
    int load_seq_from_file();

private:
    TM1Sequences<TRNAString> mMiRNASeqs;
    TM1Sequences<TRNAString> mMRNASeqs;
};

//
// TargetScan context score process core
//
template <class TRNAString, int SEEDLEN=6>
class TM1Core
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

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

public:
    // Define methods
    TM1Core(TCharSet const& pMiRNAIds, TRNASet const& pMiRNASeqs, TCharSet const& pMRNAIds, TRNASet const& pMRNASeqs,
            TIndexQGram& pRNAIdx, TFinder& pFinder):
                mExecSearchSeedSites(true), mExecGetRawFeat(true), mExecSortSites(true), mExecGetMRNAFeat(true),
                mExecRNAScore(true), mExecSumScores(true), mOutputSitePos(true), mOutputScore(true), mOutputAlign(true),
                mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
                mSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

    // Method prototypes
    void init_from_args(TM1CSOptions& opts);
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

    TM1SeedSites<TRNAString> mSeedSites;
    TM1RawFeatures<TRNAString> mRawFeatures;
    TM1SortedSitePos<TRNAString> mSortedSites;
    TM1MRNAFeatures<TRNAString> mMRNAFeatures;
    TM1MRNAInputVector<TRNAString> mMRNAInput;
    TM1ClassifiedScores<TRNAString> mScores;

private:
    int write_site_positions(seqan::CharString const &pMiRNAId);
    int write_scores(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tm1p

#endif /* TM1_CORE_HPP_ */
