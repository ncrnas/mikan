#include <seqan/seq_io.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_with_sites.hpp"  // MKRMAWithSites

using namespace seqan;

namespace mikan {

//
// MKSiteCluster methods
//
void MKRMAWithSites::clear_maps() {
    clear(mUniqRNAPosSet);
    clear(mRNASitePosMap);
    clear(mEffectiveRNAs);
}

void MKRMAWithSites::create_mrna_site_map(mikan::MKSeedSites &pSeedSites) {

    TMRNAPosSet const &mRNAPos = pSeedSites.get_mrna_pos();
    TSitePosSet const &sitePos = pSeedSites.get_site_pos();

    TSet uniqSet;
    TPosMap siteMap;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }
        uniqSet.insert((unsigned) mRNAPos[i]);
        siteMap.insert(TPosPair((unsigned) mRNAPos[i], i));
    }

    resize(mUniqRNAPosSet, uniqSet.size());
    resize(mRNASitePosMap, uniqSet.size());
    resize(mEffectiveRNAs, uniqSet.size(), true);

    unsigned idx = 0;
    TPosMap sortedPos;
    for (TItSet itSet = uniqSet.begin(); itSet != uniqSet.end(); ++itSet) {

        mUniqRNAPosSet[idx] = *itSet;

        TItMapPair itPair = siteMap.equal_range(*itSet);
        sortedPos.clear();
        unsigned count = 0;
        for (TItMap itMap = itPair.first; itMap != itPair.second; ++itMap) {
            sortedPos.insert(TPosPair(static_cast<unsigned>(sitePos[(*itMap).second]),
                                      (*itMap).second));
            ++count;
        }

        resize(mRNASitePosMap[idx], count);
        int n = 0;
        for (TItMap itMap = sortedPos.begin(); itMap != sortedPos.end(); ++itMap) {
            mRNASitePosMap[idx][n] = (*itMap).second;
            ++n;
        }

        ++idx;
    }

}

} // namespace mikan
