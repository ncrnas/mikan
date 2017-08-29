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

void MKERMAWithSites::create_pos_map() {
    resize(mUniqRNAPosSet, mUniqSetTemp.size());
    resize(mEffectiveRNAs, mUniqSetTemp.size(), true);

    unsigned idx = 0;
    for (TItSet itSet = mUniqSetTemp.begin(); itSet != mUniqSetTemp.end(); ++itSet) {

        mUniqRNAPosSet[idx] = *itSet;
        mRNAPosMap[*itSet] = idx;

        ++idx;
    }
}

} // namespace mkens

