#ifndef MKE_RNA_SCORE_HPP_
#define MKE_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_score.hpp"       // MKRNAScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_rna_score.hpp"       // MKRNAScores
#include "mke_site_score.hpp"     // MKESiteScores

namespace mkens {

//
// Total context scores
//
class MKERNAScores : public mikan::MKRNAScores {
public:
    // Define methods
    explicit MKERNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {}

};

} // namespace mkens

#endif /* MKE_RNA_SCORE_HPP_ */
