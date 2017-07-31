#ifndef PITA_SCORE_HPP_
#define PITA_SCORE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "pita_option.hpp"        // PITAOptions
#include "pita_seed_site.hpp"     // PITASeedSites
#include "vr16_ddg_core.hpp"      // VR16DDGWorkSpace

namespace ptddg {

//
// Store alignments
//
class PITAAlign {
public:
    // Constant values
    static const unsigned INDEXED_SEQ_LEN = 6;

    // Define variables
    seqan::String<bool> mEffectiveSites;
    seqan::StringSet<seqan::CharString> mAlignMRNA;
    seqan::StringSet<seqan::CharString> mAlignBars;
    seqan::StringSet<seqan::CharString> mAlignMiRNA;

public:
    // Define methods
    PITAAlign(vr16::VR16DDGWorkSpace &pVRws) : mVRws(pVRws) {}

    // Method prototype
    void clear_align();

    void resize_align(unsigned pSize);

    void create_align(int pId, mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr const &pMRNASeq,
                      seqan::CharString const &pSeedType, unsigned pSitePos, int pMismatchPos);

private:
    vr16::VR16DDGWorkSpace &mVRws;
};

//
// Store dGduplex scores
//
class PITADGDuplexScores {
public:

    // Constant values
    static const unsigned TARGET_SEQ_LEN = 50;
    static const unsigned INDEXED_SEQ_LEN = 6;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    PITADGDuplexScores(vr16::VR16DDGWorkSpace &pVRws, PITAAlign &pAlign) :
            mVRws(pVRws), mAlign(pAlign) {}

    // Method prototype
    void clear_scores();

    int calc_scores(mikan::TRNAStr const &miRNASeq, mikan::TRNASet const &pMRNASeqs, mikan::MKSeedSites &pSeedSites);

    void print_input(seqan::CharString const &pSeedType, std::string &pInputMiRNASeq, std::string &pInputMRNASeq,
                     std::vector<int> &pInputMatchSeq);

private:
    vr16::VR16DDGWorkSpace &mVRws;
    PITAAlign &mAlign;

private:
    void create_input_mirna_seq(mikan::TRNAStr const &pMiRNASeq, std::string &pInputMiRNASeq);

    void create_input_mrna_seq(mikan::TRNAStr const &pMRNASeq, int pStart, int pEnd, std::string &pInputMRNASeq);

    void create_input_matched_seq(seqan::CharString const &pSeedType, int pMismatchPos,
                                  std::vector<int> &pInputMatchSeq);
};

//
// Store dGopen scores
//
class PITADGOpenScores {
public:
    // Constant values
    static const unsigned TARGET_SEQ_LEN = 50;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const unsigned DDG_AREA = 70;
    static const unsigned DDG_OPEN = 25;

    // Define variables
    seqan::String<bool> mEffectiveSites;
    int mFlankUp;
    int mFlankDown;

    // Define methods
    PITADGOpenScores(vr16::VR16DDGWorkSpace &pVRws) : mFlankUp(0), mFlankDown(0), mVRws(pVRws) {}

    // Method prototypes
    void clear_scores();

    int calc_scores(mikan::TRNASet const &pMRNASeqs, mikan::MKSeedSites &pSeedSites);

    void print_input(std::string &pInputMRNASeq);

private:
    // Define variable
    vr16::VR16DDGWorkSpace &mVRws;

    // Method prototype
    void create_input_mrna_seq(mikan::TRNAStr const &pMRNASeq, int pStart, int pEnd, std::string &pInputMRNASeq);

};

//
// Store ddG scores
//
class PITASiteScores : public mikan::MKSiteScores {
public:
    // Define methods
    PITASiteScores(mikan::MKOptions const &opts) :
            MKSiteScores(opts),
            mDGDuplexScores(mVRws, mAlign),
            mDGOpenScores(mVRws),
            mAlign(mVRws) {}

    float const &get_score(int i) const { return mDDGScores[i]; }

    double get_dgall(int pRetIdx) { return mVRws.get_dgall(pRetIdx); }

    double get_dg5(int pRetIdx) { return mVRws.get_dg5(pRetIdx); }

    double get_dg3(int pRetIdx) { return mVRws.get_dg3(pRetIdx); }

    double get_dg0(int pRetIdx) { return mVRws.get_dg0(pRetIdx); }

    double get_dg1(int pRetIdx) { return mVRws.get_dg1(pRetIdx); }

    double get_dgsum(int pRetIdx) { return mVRws.get_dgsum(pRetIdx); }

    void set_backtrack(bool pBT) { mVRws.set_duplex_backtrack(pBT); }

    // Method prototype
    void init_from_args();

    void clear_scores();

    int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites);

    void print_alignment(int pIdx);

private:
    // Define variables
    seqan::String<float> mDDGScores;

    vr16::VR16DDGWorkSpace mVRws;
    PITADGDuplexScores mDGDuplexScores;
    PITADGOpenScores mDGOpenScores;
    PITAAlign mAlign;

};

//
// Total ddG scores
//
class PITATotalScores {
public:
    // Constant values
    const double MIN_EXP_DIFF;

    typedef std::set<unsigned>::iterator TItSet;

public:
    // Define methods
    PITATotalScores() : MIN_EXP_DIFF(-100.0) {}

    const seqan::String<float> &get_scores() { return mTotalScores; }

    const seqan::String<int> &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<int> &get_site_num() { return mSiteNum; }

    // Method prototype
    void clear_scores();

    int calc_scores(PITASeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, PITASiteScores &pDDGScores);

private:
    seqan::String<float> mTotalScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace ptddg

#endif /* PITA_SCORE_HPP_ */
