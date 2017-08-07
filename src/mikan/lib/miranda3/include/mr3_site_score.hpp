#ifndef MR3_SITE_SCORE_HPP_
#define MR3_SITE_SCORE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mr3_option.hpp"         // MR3Options
#include "mr3_align.hpp"          // MR3Align
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "vr16_fold_core.hpp"     // VR16FoldWorkSpace

namespace mr3as {

//
// Store alignment scores
//
class MR3AlignScores {
public:
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
    MR3AlignScores(vr16::VR16FoldWorkSpace &pVRws, MR3Align &pAlign) :
            mVRws(pVRws), mAlign(pAlign), mMinAlignScore(140.0) {}

    void set_min_score(float pScore) { mMinAlignScore = pScore; }

    // Method prototype
    void clear_scores();

    int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites);

    void print_input(mikan::TRNAStr &pInputMiRNASeq, mikan::TRNAStr &pInputMRNASeq);

    float get_score(int posIdx) { return mAlignScores[posIdx]; }

private:
    vr16::VR16FoldWorkSpace &mVRws;
    MR3Align &mAlign;
    seqan::String<float> mAlignScores;
    float mMinAlignScore;

private:
    void create_input_mirna_seq(mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr &pInputMiRNASeq,
                                mikan::TRNAStr &pIMiRNASeedSeq,
                                seqan::Rna5String &pIMiRNA3pSeq);

    void create_input_mrna_seq(mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr const &pMRNASeq, int pStart, int pEnd,
                               const seqan::CharString &pSeedType, mikan::TRNAStr &pInputMRNASeq,
                               mikan::TRNAStr &pIMRNASeedSeq,
                               seqan::Rna5String &pIMRNA3pSeq, bool &pNoMRNA1);

};

//
// Store energy scores
//
class MR3EnergyScores {
public:
    // Constant values
    static const unsigned TARGET_SEQ_LEN = 40;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const unsigned LINKER_LEN = 7;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MR3EnergyScores(vr16::VR16FoldWorkSpace &pVRws, MR3Align &pAlign) :
            mVRws(pVRws), mAlign(pAlign), mMaxEnergy(1.0) {}

    void set_max_score(float pScore) { mMaxEnergy = pScore; }

    float get_score(int posIdx) { return mEnScores[posIdx]; }

    // Method prototype
    void clear_scores();

    int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites);

    void print_input(std::string &pInputMRNASeq);

private:
    vr16::VR16FoldWorkSpace &mVRws;
    MR3Align &mAlign;
    seqan::String<float> mEnScores;
    float mMaxEnergy;

private:
    void create_input_seq(int pIdx, mikan::TRNAStr const &pMiRNASeq, std::string &pInputMRNASeq);
};

//
// Store miRanda scores
//
class MR3SiteScores : public mikan::MKSiteScores {
public:
    // Constant values
    static const unsigned RNAFOLD_MAX_INPUTLEN = 60;

    // Define methods
    explicit MR3SiteScores(mikan::MKOptions const &opts) :
            MKSiteScores(opts),
            mVRws(30.0),
            mAlign(),
            mAlignScores(mVRws, mAlign),
            mEnergyScores(mVRws, mAlign) {
        init_rnafold();
    }

    virtual float get_score(int pIdx) { return get_align_score(pIdx); }

    float get_align_score(int posIdx) { return mAlignScores.get_score(posIdx); }

    float get_energy_score(int posIdx) { return mEnergyScores.get_score(posIdx); }

    void set_backtrack(bool pBT) { mVRws.set_backtrack(pBT); }

    void set_min_align_score(float pScore) { mAlignScores.set_min_score(pScore); }

    void set_max_energy(float pScore) { mEnergyScores.set_max_score(pScore); }

    // Method prototype
    virtual void clear_scores();

    virtual int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites &pRNAWithSites);


    void print_alignment(int pIdx);

    void init_rnafold();

private:
    vr16::VR16FoldWorkSpace mVRws;
    MR3Align mAlign;
    MR3AlignScores mAlignScores;
    MR3EnergyScores mEnergyScores;

};

} // namespace mr3as

#endif /* MR3_SITE_SCORE_HPP_ */
