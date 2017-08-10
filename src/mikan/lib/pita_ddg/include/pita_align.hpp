#ifndef PITA_ALIGN_HPP_
#define PITA_ALIGN_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "pita_option.hpp"        // PITAOptions
#include "pita_seed_site.hpp"     // PITASeedSites
#include "vr16_ddg_core.hpp"      // VR16DDGWorkSpace

namespace ptddg {

//
// Store alignments
//
class PITAAlign {
public:
    // Define variables
    seqan::String<bool> mEffectiveSites;
    seqan::StringSet<seqan::CharString> mAlignMRNA;
    seqan::StringSet<seqan::CharString> mAlignBars;
    seqan::StringSet<seqan::CharString> mAlignMiRNA;

public:
    // Define methods
    explicit PITAAlign(vr16::VR16DDGWorkSpace &pVRws) : mVRws(pVRws) {}

    // Method prototype
    void clear_align();

    void resize_align(unsigned pSize);

    void create_align(int pId, mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr const &pMRNASeq,
                      seqan::CharString const &pSeedType, unsigned pSitePos, int pMismatchPos);

private:
    vr16::VR16DDGWorkSpace &mVRws;
};

} // namespace ptddg

#endif /* PITA_ALIGN_HPP_ */
