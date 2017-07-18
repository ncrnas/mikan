#ifndef MR3_SEED_SITE_HPP_
#define MR3_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"        // MKSeedSeqs
#include "mk_seed_site.hpp"       // MKSeedSites

namespace mr3as {

//
// Generate miRNA seeds
//
class MR3SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    MR3SeedSeqs() : MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);
};

//
// miRNA seed sites
//
class MR3SeedSites : public mikan::MKSeedSites {
public:
    // Constant values
    static const bool FORCE_LAST_MATCH = false;

public:
    // Define methods
    MR3SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {
        mMinToCDS = 1;
        mMinToEnd = 6;
    }

private:
    bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                           mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                           seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

    void set_mx_matches(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, int pMx,
                        bool &pMatchMx, bool &pGutMx, bool &pGumMx);

    void set_stringent_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 bool pMatchMx8, bool pMatchMx9, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_single_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                 bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_multiple_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                   bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void set_gu_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                   bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_bt_seed_type(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, unsigned pMisMatchPos,
                          seqan::CharString &pNewSeedType);

    void set_6mer_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                            bool pMatchMx8, bool pMatchMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType);

};

} // namespace mr3as

#endif /* MR3_SEED_SITE_HPP_ */
