#ifndef MKE_SEED_SITE_HPP_
#define MKE_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs
#include "mk_seed_site.hpp"      // MKSeedSites
#include "mk_option.hpp"         // MKOptions
#include "mke_option.hpp"        // MKEOptions

namespace mkens {

//
// miRNA seed sites
//
class MKESeedSites : public mikan::MKSeedSites {
public:
    // Define methods
    MKESeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {

        resize(mSeedTypeList, mkens::TOOL_NUM);

    }

    unsigned get_idx_from_pos(unsigned pMRNAPos, unsigned pSitePos) {
        TPosPair pair = std::make_pair(pMRNAPos, pSitePos);
        return mPosPairMap[pair];
    };

    // Method prototypes
    virtual void clear_pos();

    void add_to_set(mikan::MKSeedSites &pSeedSites, unsigned pToolIdx, seqan::CharString &pPrefix);
    void add_seed_types(mikan::MKSeedSites &pSeedSites, unsigned pToolIdx, seqan::CharString &pPrefix);
    void set_default_seed_type(unsigned pIdx, seqan::CharString &pPrefix);
    void create_pos_map();
    void combine_seed_types();

private:
    // Define types
    typedef std::pair<unsigned, unsigned> TPosPair;
    typedef std::set<TPosPair> TSet;
    typedef std::set<TPosPair>::iterator TItSet;
    typedef std::map<TPosPair, unsigned> TPosMap;

    // Define variables
    TSet mUniqSet;
    TPosMap mPosPairMap;
    seqan::StringSet<seqan::StringSet<seqan::CharString> > mSeedTypeList;

};

} // namespace mkens

#endif /* MKE_SEED_SITE_HPP_ */
