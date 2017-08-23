#ifndef MKE_CORE_HPP_
#define MKE_CORE_HPP_

#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_core_base.hpp"       // MKCoreBase
#include "mke_option.hpp"         // MKEOptions
#include "mke_site_score.hpp"     // MKESiteScores
#include "mke_seed_site.hpp"      // MKESeedSites
#include "mke_rna_score.hpp"      // MKERNAScores
#include "mr3_core.hpp"           // MR3Core
#include "pita_core.hpp"          // PITACore
#include "rh2_core.hpp"           // RH2Core
#include "tm1_core.hpp"           // TM1Core
#include "ts5_core.hpp"           // TS5Core
#include "tssvm_core.hpp"         // TSSVMCore

namespace mkens {

//
// mikan ensemble score process core
//

class MKECore : public mikan::MKCoreBase  {
public:
    // Define methods
    MKECore(mkens::MKEOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
            mSiteScores(pOpts),
            mRNAScores(pOpts),
            mMR3Core(pOpts.mMR3Opts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mPITACore(pOpts.mPITAOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mRH2Core(pOpts.mRH2Opts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mTM1Core(pOpts.mTM1Opts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mTS5Core(pOpts.mTS5Opts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mTSSVMCore(pOpts.mTSSVMOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {}

    virtual int find_seed_sites(unsigned pIdx);
    virtual int calc_site_scores(unsigned pIdx);
    virtual int calc_rna_scores(unsigned);
    virtual int output_results(unsigned pIdx);
    virtual void clear_all();

    virtual mikan::MKSeedSites & get_seed_sites() { return mSeedSites; }
    virtual mikan::MKSiteScores & get_site_scores() { return mSiteScores; };
    virtual mikan::MKRNAScores & get_rna_scores() { return mRNAScores; };

private:
    MKESeedSites mSeedSites;
    MKESiteScores mSiteScores;
    MKERNAScores mRNAScores;

    mr3as::MR3Core mMR3Core;
    ptddg::PITACore mPITACore;
    rh2mfe::RH2Core mRH2Core;
    tm1p::TM1Core mTM1Core;
    ts5cs::TS5Core mTS5Core;
    tssvm::TSSVMCore mTSSVMCore;

    mikan::MKCoreBase& get_tool_core(unsigned pIdx) {
        switch (pIdx) {
            case 0:
                return mMR3Core;
            case 1:
                return mPITACore;
            case 2:
                return mRH2Core;
            case 3:
                return mTM1Core;
            case 4:
                return mTS5Core;
            case 5:
                return mTSSVMCore;
            default:
                return *this;
        }
    }

    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);

    int combine_site_pos(unsigned pIdx);
    int combine_seed_types();

};

} // namespace mkens

#endif /* MKE_CORE_HPP_ */
