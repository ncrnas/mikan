#ifndef PITA_SEED_SITE_HPP_
#define PITA_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"        // MKSeedSeqs
#include "mk_seed_site.hpp"       // MKSeedSites

namespace ptddg {

//
// Generate miRNA seeds
//
class PITASeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    PITASeedSeqs() : MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);

};

//
// miRNA seed sites
//
class PITASeedSites : public mikan::MKSeedSites {
public:
    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 22;
    static const unsigned MIN_DIST_UTR_END = 1;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const bool FORCE_LAST_MATCH = true;

    // Define methods
    PITASeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

private:
    virtual bool check_position(unsigned pMRNAPos, unsigned pSitePos);

    virtual void set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                                   mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                                   seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

    void set_mx_matches(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, int pMx,
                        bool &pMatchMx, bool &pGutMx, bool &pGumMx);

    void set_stringent_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 bool pMatchMx1, bool pMatchMx2, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_single_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                 bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_multiple_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                   bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void set_gu_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                   bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_6mer_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                            bool pMatchMx1, bool pMatchMx2, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType);

};

} // namespace ptddg

#endif /* PITA_SEED_SITE_HPP_ */
