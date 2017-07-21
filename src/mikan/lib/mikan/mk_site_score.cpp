#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"     // MKSiteScores

using namespace seqan;

namespace mikan {

//
// MKSiteScores methods
//
void MKSiteScores::clear_scores() {
    clear(mEffectiveSites);
    clear(mSiteScores);
}

void MKSiteScores::init_from_args(MKOptions &) {

}

int MKSiteScores::calc_scores(
        mikan::TRNAStr const &,
        mikan::TRNASet const &,
        mikan::MKSeedSites &pSeedSites) {

    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    resize(mSiteScores, length(pSeedSites.mEffectiveSites));

    for (unsigned i = 0; i < length(mEffectiveSites); i++) {
        if (!pSeedSites.mEffectiveSites[i]) {
            mSiteScores[i] = 0;
            mEffectiveSites[i] = false;
            continue;
        }

        mSiteScores[i] = 1;
        mEffectiveSites[i] = true;
    }

    return 0;
}

} // namespace mikan
