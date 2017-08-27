#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mke_seed_site.hpp"     // MKESeedSeqs, MKESeedSites

using namespace seqan;

namespace mkens {

//
// MKESeedSites methods
//
void MKESeedSites::clear_pos() {
    mikan::MKSeedSites::clear_pos();

    clear(mUniqSet);
    clear(mPosPairMap);

    for (unsigned i = 0; i < length(mSeedTypeList); i++) {
        clear(mSeedTypeList[i]);
    }
}

void MKESeedSites::add_to_set(
        mikan::MKSeedSites &pSeedSites,
        unsigned,
        seqan::CharString &) {
    const mikan::TMRNAPosSet &RNAPos = pSeedSites.get_mrna_pos();
    const mikan::TMRNAPosSet &S1Pos = pSeedSites.get_site_pos_s1();

    for (unsigned i = 0; i < length(pSeedSites.mEffectiveSites); i++) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }
        TPosPair pair = std::make_pair(RNAPos[i], S1Pos[i]);
        mUniqSet.insert(pair);
    }
}

void MKESeedSites::create_pos_map() {
    resize(mEffectiveSites, mUniqSet.size(), true);
    resize(mMRNAPos, mUniqSet.size());
    resize(mSitePos, mUniqSet.size());
    resize(mSeedTypes, mUniqSet.size());
    resize(mMisMatchPos, mUniqSet.size(), 0);

    unsigned idx = 0;
    for (TItSet itSet = mUniqSet.begin(); itSet != mUniqSet.end(); ++itSet) {

        mMRNAPos[idx] = (*itSet).first;
        mSitePos[idx] = (*itSet).second;
        mSeedTypes[idx] = "";
        mPosPairMap[*itSet] = idx;

        ++idx;
    }
}

void MKESeedSites::add_seed_types(
        mikan::MKSeedSites &pSeedSites,
        unsigned pToolIdx,
        seqan::CharString &pPrefix) {

    const mikan::TMRNAPosSet &RNAPos = pSeedSites.get_mrna_pos();
    const mikan::TMRNAPosSet &S1Pos = pSeedSites.get_site_pos_s1();
    const seqan::StringSet<seqan::CharString> &seedTypes = pSeedSites.get_seed_types();

    resize(mSeedTypeList[pToolIdx], length(mEffectiveSites));
    set_default_seed_type(pToolIdx, pPrefix);

    seqan::CharString prefix;
    append(prefix, pPrefix);
    append(prefix, ":");

    for (unsigned i = 0; i < length(pSeedSites.mEffectiveSites); i++) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }

        unsigned idx = get_idx_from_pos(RNAPos[i], S1Pos[i]);
        seqan::CharString seedType;
        append(seedType, pPrefix);
        append(seedType, ":");
        append(seedType, seedTypes[i]);

        mSeedTypeList[pToolIdx][idx] = seedType;
    }
}

void MKESeedSites::set_default_seed_type(unsigned pIdx, seqan::CharString &pPrefix) {
    seqan::CharString defSeedType;
    append(defSeedType, pPrefix);
    append(defSeedType, ":NA");

    for (unsigned i = 0; i < length(mEffectiveSites); i++) {
        append(mSeedTypeList[pIdx][i], defSeedType);
    }

}

void MKESeedSites::combine_seed_types() {
    for (unsigned i = 0; i < length(mEffectiveSites); i++) {
        if (!mEffectiveSites[i]) {
            continue;
        }

        seqan::CharString seedType = "";
        for (unsigned j = 0; j < length(mSeedTypeList); j++) {
            append(seedType, mSeedTypeList[j][i]);
            append(seedType, ",");
        }

        mSeedTypes[i] = seedType;
    }
}

} // namespace mkens

