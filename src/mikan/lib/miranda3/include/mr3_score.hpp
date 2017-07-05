#ifndef MR3_SCORE_HPP_
#define MR3_SCORE_HPP_

#include <mr3_align.hpp>          // MR3SeedSites
#include <mr3_seed_site.hpp>      // MR3Align
#include <vr16_fold_core.hpp>     // VR16FoldWorkSpace
#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>

namespace mr3as {

//
// Store alignment scores
//
template<class TRNAString>
class MR3AlignScores {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Constant values
    static const unsigned TARGET_SEQ_LEN = 40;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const unsigned SEED_REGION_LEN = 8;
    static const int MATCH_SCORE = 5;
    static const int EXTENT_SCORE = -4;
    static const int OFFSET3P = 2;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MR3AlignScores(vr16::VR16FoldWorkSpace &pVRws, MR3Align<TRNAString> &pAlign) :
            mVRws(pVRws), mAlign(pAlign), mMinAlignScore(140.0) {}

    void set_min_score(float pScore) { mMinAlignScore = pScore; }

    // Method prototype
    void clear_scores();

    int calc_scores(MR3SeedSites<TRNAString> &pSeedSites, TRNAString const &miRNASeq, TRNASet const &pMRNASeqs);

    void print_input(TRNAString &pInputMiRNASeq, TRNAString &pInputMRNASeq);

    float get_score(int posIdx) { return mAlignScores[posIdx]; }

private:
    vr16::VR16FoldWorkSpace &mVRws;
    MR3Align<TRNAString> &mAlign;
    seqan::String<float> mAlignScores;
    float mMinAlignScore;

private:
    void create_input_mirna_seq(TRNAString const &pMiRNASeq, TRNAString &pInputMiRNASeq, TRNAString &pIMiRNASeedSeq,
                                seqan::Rna5String &pIMiRNA3pSeq);

    void create_input_mrna_seq(TRNAString const &pMiRNASeq, TRNAString const &pMRNASeq, int pStart, int pEnd,
                               const seqan::CharString &pSeedType, TRNAString &pInputMRNASeq,
                               TRNAString &pIMRNASeedSeq,
                               seqan::Rna5String &pIMRNA3pSeq, bool &pNoMRNA1);

};

//
// Store energy scores
//
template<class TRNAString>
class MR3EnergyScores {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Constant values
    static const unsigned TARGET_SEQ_LEN = 40;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const unsigned LINKER_LEN = 7;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MR3EnergyScores(vr16::VR16FoldWorkSpace &pVRws, MR3Align<TRNAString> &pAlign) :
            mVRws(pVRws), mAlign(pAlign), mMaxEnergy(1.0) {}

    void set_max_score(float pScore) { mMaxEnergy = pScore; }

    float get_score(int posIdx) { return mEnScores[posIdx]; }

    // Method prototype
    void clear_scores();

    int calc_scores(MR3SeedSites<TRNAString> &pSeedSites, TRNAString const &miRNASeq, TRNASet const &pMRNASeqs);

    void print_input(std::string &pInputMRNASeq);

private:
    vr16::VR16FoldWorkSpace &mVRws;
    MR3Align<TRNAString> &mAlign;
    seqan::String<float> mEnScores;
    float mMaxEnergy;

private:
    void create_input_seq(int pIdx, TRNAString const &pMiRNASeq, std::string &pInputMRNASeq);
};

//
// Store miRanda scores
//
template<class TRNAString>
class MR3SiteScores {
public:
    // Constant values
    static const unsigned RNAFOLD_MAX_INPUTLEN = 60;

    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MR3SiteScores() : mVRws(30.0), mAlign(), mAlignScores(mVRws, mAlign), mEnergyScores(mVRws, mAlign) {
        init_rnafold();
    }

    float get_align_score(int posIdx) { return mAlignScores.get_score(posIdx); }

    float get_energy_score(int posIdx) { return mEnergyScores.get_score(posIdx); }

    void set_backtrack(bool pBT) { mVRws.set_backtrack(pBT); }

    void set_min_align_score(float pScore) { mAlignScores.set_min_score(pScore); }

    void set_max_energy(float pScore) { mEnergyScores.set_max_score(pScore); }

    // Method prototype
    void clear_scores();

    int calc_scores(MR3SeedSites<TRNAString> &pSeedSites, TRNAString const &miRNASeq, TRNASet const &pMRNASeqs);

    void print_alignment(int pIdx);

    void init_rnafold();

private:
    vr16::VR16FoldWorkSpace mVRws;
    MR3Align<TRNAString> mAlign;
    MR3AlignScores<TRNAString> mAlignScores;
    MR3EnergyScores<TRNAString> mEnergyScores;

};

//
// Total scores
//
template<class TRNAString>
class MR3TotalScores {
public:
    // Constant values
    const double MIN_EXP_DIFF;

public:
    // Define methods
    MR3TotalScores() : MIN_EXP_DIFF(-100.0) {}

    const seqan::String<float> &get_align_scores() { return mTotalAlignScores; }

    const seqan::String<float> &get_energy_scores() { return mTotalEnScores; }

    const seqan::String<int> &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<int> &get_site_num() { return mSiteNum; }

    // Method prototype
    void clear_scores();

    int calc_scores(MR3SeedSites<TRNAString> &pSeedSites, MR3SiteScores<TRNAString> &pSiteScores,
                    const seqan::String<unsigned> &pSortedSites);

private:
    seqan::String<float> mTotalAlignScores;
    seqan::String<float> mTotalEnScores;
    seqan::String<float> mLogMaxAlignScores;
    seqan::String<float> mLogMaxEnScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace mr3as

#endif /* MR3_SCORE_HPP_ */
