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

void MKRMAWithSites::create_temp_map(mikan::MKSeedSites &pSeedSites) {
    TMRNAPosSet const &mRNAPos = pSeedSites.get_mrna_pos();

    mUniqSetTemp.clear();
    mSiteMapTemp.clear();
    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }
        mUniqSetTemp.insert((unsigned) mRNAPos[i]);
        mSiteMapTemp.insert(TPosPair((unsigned) mRNAPos[i], i));
    }
}

void MKRMAWithSites::create_mrna_site_map(
        mikan::MKSeedSites &pSeedSites,
        MKSiteScores &pSiteScores) {

    TSitePosSet const &sitePos = pSeedSites.get_site_pos();

    create_temp_map(pSeedSites);

    resize(mUniqRNAPosSet, mUniqSetTemp.size());
    resize(mRNASitePosMap, mUniqSetTemp.size());
    resize(mEffectiveRNAs, mUniqSetTemp.size(), true);

    unsigned idx = 0;
    TPosMap sortedPos;
    for (TItSet itSet = mUniqSetTemp.begin(); itSet != mUniqSetTemp.end(); ++itSet) {

        mUniqRNAPosSet[idx] = *itSet;

        TItMapPair itPair = mSiteMapTemp.equal_range(*itSet);
        sortedPos.clear();
        unsigned count = 0;
        for (TItMap itMap = itPair.first; itMap != itPair.second; ++itMap) {
            if (mSortVtype == "score") {

                sortedPos.insert(TPosPair(pSiteScores.get_score((*itMap).second),
                                          (*itMap).second));
            } else if (mSortVtype == "wide") {
                sortedPos.insert(TPosPair(static_cast<float>(pSiteScores.get_wide_site_start((*itMap).second)),
                                          (*itMap).second));
            } else {
                sortedPos.insert(TPosPair(static_cast<float>(sitePos[(*itMap).second]),
                                          (*itMap).second));
            }

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
