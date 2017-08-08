#ifndef PITA_SEED_SITE_HPP_
#define PITA_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"        // MKSeedSeqs
#include "mk_seed_site.hpp"       // MKSeedSites
#include "mk_option.hpp"          // MKOptions

namespace ptddg {

//
// Generate miRNA seeds
//
class PITASeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    PITASeedSeqs(mikan::MKOptions const &opts) : MKSeedSeqs(opts) {
        init_from_args();
        set_flags();
    }

    // Method prototypes
    void init_from_args();
    void set_flags();
};

//
// miRNA seed sites
//
class PITASeedSites : public mikan::MKSeedSites {
public:
    // Constant values
    static const bool FORCE_LAST_MATCH = true;

    // Define methods
    PITASeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {
        
        mMinToCDS = 16;
        mMinToEnd = 7;
    }

private:
    virtual bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                           mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                           seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

    void set_stringent_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 bool pMatchM8, bool pMatchM9, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_single_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 int pM1, int pM2, bool pMatchM8, bool pMatchM9, bool pGutM8, bool pGutM9,
                                 bool pGumM8, bool pGumM9, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_multiple_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchM8, bool pMatchM9, bool pGutM8, bool pGutM9,
                                   bool pGumM8, bool pGumM9, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                int pM1, int pM2, bool pMatchM8, bool pMatchM9, bool pGutM8, bool pGutM9,
                                bool pGumM8, bool pGumM9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void set_gu_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchM8, bool pMatchM9, bool pGutM8, bool pGutM9,
                                   bool pGumM8, bool pGumM9, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_6mer_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                            bool pMatchM8, bool pMatchM9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType);

};

} // namespace ptddg

#endif /* PITA_SEED_SITE_HPP_ */
