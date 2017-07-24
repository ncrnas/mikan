#include <seqan/seq_io.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_with_sites.hpp"  // MKRMAWithSites

using namespace seqan;

namespace mikan {

//
// MKSiteCluster methods
//
void MKRMAWithSites::clear_sets() {
    mEffectiveSiteCount = 0;
    mUniqRNASet.clear();
    mRNASiteMap.clear();
    clear(mSortedMRNAPos);
}

void MKRMAWithSites::cluster_sites(mikan::MKSeedSites &pSeedSites) {

    String<unsigned> const &mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }
        mUniqRNASet.insert((unsigned) mRNAPos[i]);
        mRNASiteMap.insert(TPosPair((unsigned) mRNAPos[i], i));
        ++mEffectiveSiteCount;
    }
}

void MKRMAWithSites::sort_mrna_pos(MKSeedSites &pSeedSites) {
    clear_sets();
    cluster_sites(pSeedSites);

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();

    resize(mSortedMRNAPos, mUniqRNASet.size());

    TPosMap sortedPos;
    unsigned idx = 0;
    for (TItSet itSet = mUniqRNASet.begin(); itSet != mUniqRNASet.end(); ++itSet) {

        TItMapPair itPair = mRNASiteMap.equal_range(*itSet);
        sortedPos.clear();
        unsigned count = 0;

        for (TItMap itMap = itPair.first; itMap != itPair.second; ++itMap) {
            sortedPos.insert(TPosPair(static_cast<unsigned>(sitePos[(*itMap).second]),
                                      (*itMap).second));
            ++count;
        }

        resize(mSortedMRNAPos[idx], count);
        int n = 0;
        for (TItMap itMap = sortedPos.begin(); itMap != sortedPos.end(); ++itMap) {
            mSortedMRNAPos[idx][n] = (*itMap).second;
            ++n;
        }

        ++idx;
    }
}

} // namespace mikan
