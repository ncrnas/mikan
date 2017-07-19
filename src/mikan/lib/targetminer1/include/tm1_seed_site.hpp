#ifndef TM1_SEED_SITE_HPP_
#define TM1_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs
#include "mk_seed_site.hpp"      // MKSeedSites

namespace tm1p {

//
// Generate miRNA seeds
//
class TM1SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    TM1SeedSeqs() : MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);

protected:
    virtual int create_other_seed_seqs(mikan::TRNAStr &pSeedSeq);
};

//
// miRNA seed sites
//
class TM1SeedSites : public mikan::MKSeedSites {
public:
    // Define methods
    TM1SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

    bool is_m8_match_gu(int i) { return (mM8Match[i] || mM8GU[i]); }

    bool is_m8_match(int i) { return (mM8Match[i]); }

    // Method prototypes

    int get_seed_len(int pIdx);

    int get_seed_start_pos(int pIdx);

    int get_seed_start_pos2(int pIdx);

    int get_seed_end_pos(int pIdx);

    int get_seed_end_pos2(int pIdx);

    int get_length_to_cds(int pIdx);

    void clear_pos();

private:
    seqan::String<bool> mM8Match;
    seqan::String<bool> mM8GU;
    seqan::String<bool> mM1A;
    seqan::String<bool> mM1Match;
    seqan::String<bool> mM1GU;
    seqan::String<unsigned> mMRNASeqLen;

private:
    bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                           mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                           seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

    void get_match_count(unsigned pSitePos, unsigned pMRNAPos, mikan::TRNAStr const &pMiRNASeq,
                         int pMx1, int pMx2, int &pMatchCount, int &pGUCount);
};


} // namespace tm1p

#endif /* TM1_SEED_SITE_HPP_ */
