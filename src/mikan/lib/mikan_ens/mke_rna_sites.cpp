#include <seqan/seq_io.h>
#include "mke_rna_sites.hpp"     // MKERMAWithSites

using namespace seqan;

namespace mkens {

//
// MKERMAWithSites methods
//
void MKERMAWithSites::clear_maps() {
    mikan::MKRMAWithSites::clear_maps();

    mRNAPosMap.clear();
}

void MKERMAWithSites::add_to_set(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &pSiteScores) {

    const mikan::TMRNAPosSet &RNAPos = pSeedSites.get_mrna_pos();
    mUniqSetTemp.clear();
    mSiteMapTemp.clear();

    for (unsigned i = 0; i < length(pSeedSites.mEffectiveSites); i++) {
        if (!pSeedSites.mEffectiveSites[i] || !pSiteScores.mEffectiveSites[i]) {
            continue;
        }
        mUniqSetTemp.insert((unsigned) RNAPos[i]);
        mSiteMapTemp.insert(TPosPair((unsigned) RNAPos[i], i));
    }
}

void MKERMAWithSites::create_pos_map(mikan::MKSeedSites &pSeedSites) {
    mikan::TSitePosSet const &sitePos = pSeedSites.get_site_pos();

    resize(mUniqRNAPosSet, mUniqSetTemp.size());
    resize(mRNASitePosMap, mUniqSetTemp.size());
    resize(mEffectiveRNAs, mUniqSetTemp.size(), true);

    unsigned idx = 0;
    TPosMap sortedPos;
    for (TItSet itSet = mUniqSetTemp.begin(); itSet != mUniqSetTemp.end(); ++itSet) {

        mUniqRNAPosSet[idx] = *itSet;
        mRNAPosMap[*itSet] = idx;

        TItMapPair itPair = mSiteMapTemp.equal_range(*itSet);
        sortedPos.clear();
        unsigned count = 0;
        for (TItMap itMap = itPair.first; itMap != itPair.second; ++itMap) {
            sortedPos.insert(TPosPair(static_cast<float>(sitePos[(*itMap).second]),
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

} // namespace mkens

