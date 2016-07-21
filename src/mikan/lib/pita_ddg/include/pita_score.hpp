#ifndef PITA_SCORE_HPP_
#define PITA_SCORE_HPP_

#include <mikan/lib/pita_ddg/include/pita_seed_site.hpp>     // PITASeedSites
#include <mikan/lib/vienna_rna/include/vr16_ddg_core.hpp>      // VR16DDGWorkSpace
#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>

namespace ptddg{

//
// Store alignments
//
template <class TRNAString>
class PITAAlign
{
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
    PITAAlign(vr16::VR16DDGWorkSpace& pVRws): mVRws(pVRws) {}

    // Method prototype
    void clear_align();
    void resize_align(unsigned pSize);
    void create_align(int pId, TRNAString const &pMiRNASeq, TRNAString const &pMRNASeq,
            seqan::CharString const &pSeedType, unsigned pSitePos, int pMismatchPos);

private:
    vr16::VR16DDGWorkSpace& mVRws;
};

//
// Store dGduplex scores
//
template <class TRNAString>
class PITADGDuplexScores
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Constant values
    static const unsigned TARGET_SEQ_LEN = 50;
    static const unsigned INDEXED_SEQ_LEN = 6;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    PITADGDuplexScores(vr16::VR16DDGWorkSpace& pVRws, PITAAlign<TRNAString>& pAlign) :
        mVRws(pVRws), mAlign(pAlign) {}

    // Method prototype
    void clear_scores();
    int calc_scores(PITASeedSites<TRNAString> &pSeedSites, TRNAString const &miRNASeq, TRNASet const &pMRNASeqs);
    void print_input(seqan::CharString const &pSeedType, std::string &pInputMiRNASeq, std::string &pInputMRNASeq,
            std::vector<int> &pInputMatchSeq);

private:
    vr16::VR16DDGWorkSpace& mVRws;
    PITAAlign<TRNAString>& mAlign;

private:
    void create_input_mirna_seq(TRNAString const &pMiRNASeq, std::string &pInputMiRNASeq);
    void create_input_mrna_seq(TRNAString const &pMRNASeq, int pStart, int pEnd, std::string &pInputMRNASeq);
    void create_input_matched_seq(seqan::CharString const &pSeedType, int pMismatchPos,
            std::vector<int> &pInputMatchSeq);
};

//
// Store dGopen scores
//
template <class TRNAString>
class PITADGOpenScores
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Constant values
    static const unsigned TARGET_SEQ_LEN = 50;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const unsigned DDG_AREA = 70;
    static const unsigned DDG_OPEN = 25;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    PITADGOpenScores(vr16::VR16DDGWorkSpace& pVRws) : mVRws(pVRws) {}

    // Method prototype
    void clear_scores();
    int calc_scores(PITASeedSites<TRNAString> &pSeedSites, TRNASet const &pMRNASeqs, int pFlankUp, int pFlankDown);
    void print_input(std::string &pInputMRNASeq);

private:
    vr16::VR16DDGWorkSpace &mVRws;

private:
    void create_input_mrna_seq(TRNAString const &pMRNASeq, int pStart, int pEnd, std::string &pInputMRNASeq);
};

//
// Store ddG scores
//
template <class TRNAString>
class PITADDGScores
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    PITADDGScores() : mDGDuplexScores(mVRws, mAlign), mDGOpenScores(mVRws), mAlign(mVRws) {}
    float const& get_score(int i) const {return mDDGScores[i];}
    double get_dgall(int pRetIdx) {return mVRws.get_dgall(pRetIdx);}
    double get_dg5(int pRetIdx) {return mVRws.get_dg5(pRetIdx);}
    double get_dg3(int pRetIdx) {return mVRws.get_dg3(pRetIdx);}
    double get_dg0(int pRetIdx) {return mVRws.get_dg0(pRetIdx);}
    double get_dg1(int pRetIdx) {return mVRws.get_dg1(pRetIdx);}
    double get_dgsum(int pRetIdx) {return mVRws.get_dgsum(pRetIdx);}
    void set_backtrack(bool pBT){mVRws.set_duplex_backtrack(pBT);}

    // Method prototype
    void clear_scores();
    int calc_scores(PITASeedSites<TRNAString> &pSeedSites, TRNAString const &miRNASeq, TRNASet const &pMRNASeqs,
            int pFlankUp, int pFlankDown);
    void print_alignment(int pIdx);

private:
    seqan::String<float> mDDGScores;

    vr16::VR16DDGWorkSpace mVRws;
    PITADGDuplexScores<TRNAString> mDGDuplexScores;
    PITADGOpenScores<TRNAString> mDGOpenScores;
    PITAAlign<TRNAString> mAlign;

};

//
// Total ddG scores
//
template <class TRNAString>
class PITATotalScores
{
public:
    // Constant values
    const double MIN_EXP_DIFF;

public:
    // Define methods
    PITATotalScores(): MIN_EXP_DIFF(-100.0) {}
    const seqan::String<float>& get_scores(){return mTotalScores;}
    const seqan::String<int>& get_mrna_pos(){return mMRNAPos;}
    const seqan::String<int>& get_site_num(){return mSiteNum;}

    // Method prototype
    void clear_scores();
    int calc_scores(PITASeedSites<TRNAString> &pSeedSites, PITADDGScores<TRNAString> &pDDGScores,
            const seqan::String<unsigned> &pSortedSites);

private:
    seqan::String<float> mTotalScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace ptddg

#endif /* PITA_SCORE_HPP_ */
