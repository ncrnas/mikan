#ifndef MKE_RNA_SITES_HPP_
#define MKE_RNA_SITES_HPP_

#include <seqan/sequence.h>
#include "mk_rna_sites.hpp"         // MKRMAWithSites

namespace mkens {

//
// Cluster target sites
//
class MKERMAWithSites : public mikan::MKRMAWithSites {
public:
    // Define methods
    explicit MKERMAWithSites(mikan::MKOptions const &opts) : MKRMAWithSites(opts) {}

    unsigned get_idx_from_pos(unsigned pMRNAPos) { return mRNAPosMap[pMRNAPos]; }

    // Method prototypes
    virtual void clear_maps();

    void add_to_set(mikan::MKSeedSites &pSeedSites, mikan::MKSiteScores &pSiteScores);

    void create_pos_map(mikan::MKSeedSites &pSeedSites);

private:
    typedef std::map<unsigned, unsigned> TPosMapU;

    TPosMapU mRNAPosMap;

};

} // namespace mkens

#endif /* MKE_RNA_SITES_HPP_ */
