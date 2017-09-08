#ifndef TSSVM_SITE_SCORE_HPP_
#define TSSVM_SITE_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_option.hpp"          // MKOptions
#include "mk_seed_site.hpp"       // MKSeedSites
#include "tssvm_align.hpp"        // TSSVMAlign
#include "tssvm_option.hpp"       // TSSVMOptions
#include "tssvm_seed_site.hpp"    // TSSVMSeedSites, TSSVMSiteFilter
#include "tssvm_site_svm.hpp"     // TSSVMSiteInputVector

namespace tssvm {

//
// Store site level scores
//
class TSSVMSiteScores : public mikan::MKSiteScores {
public:
    // Define methods
    explicit TSSVMSiteScores(mikan::MKOptions const &opts) :
            MKSiteScores(opts),
            mSiteInput(mSiteModel) {

        resize(mScoreTypes, 1);
        mScoreTypes[0] = "svm";
    }

    virtual float get_score(int i) { return mSiteInput.get_score(i); }

    virtual float get_score(int, int pIdx) { return mSiteInput.get_score(pIdx); }

    const seqan::String<float> &get_scores() { return mSiteInput.get_scores(); }

    void write_alignment(int pIdx) { mAlignSeqs.write_alignment(pIdx); }

    // Method prototypes
    virtual void clear_scores();

    virtual int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites &pRNAWithSites);

private:
    TSSVMAlign mAlignSeqs;
    TSSVMRawFeatures mSiteFeatures;
    TSSVMSiteModel mSiteModel;
    TSSVMSiteInputVector mSiteInput;

};

} // namespace tssvm

#endif /* TSSVM_SITE_SCORE_HPP_ */
