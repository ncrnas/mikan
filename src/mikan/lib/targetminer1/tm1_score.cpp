#include <math.h>                // roundf
#include <mikan/lib/targetminer1/include/tm1_inst_template.hpp> // TRNATYPE
#include <mikan/lib/targetminer1/include/tm1_score.hpp>         // TM1ClassifiedScores
#include <seqan/sequence.h>

using namespace seqan;

namespace tm1p{

//
// TM1TotalScores methods
//
template <class TRNAString>
void TM1ClassifiedScores<TRNAString>::clear_scores()
{
    clear(mScores);
    clear(mPredictions);
    clear(mSiteNum);
}

template <class TRNAString>
int TM1ClassifiedScores<TRNAString>::calc_scores(
        const seqan::String<unsigned>& pSiteCoutns,
        const seqan::String<float>& pScores)
{
    resize(mScores, length(pSiteCoutns));
    resize(mPredictions, length(pSiteCoutns));
    resize(mSiteNum, length(pSiteCoutns));

    for (unsigned i = 0; i < length(pSiteCoutns); ++i)
    {

        mScores[i] = roundf(pScores[i] * 1000.0 ) / 1000.0;
        if (pScores[i] > 0)
        {
            mPredictions[i]= 1;
        }
        else
        {
            mPredictions[i]= -1;
        }
        mSiteNum[i] = pSiteCoutns[i];
    }

    return 0;
}

// Explicit template instantiation
template class TM1ClassifiedScores<TRNATYPE>;

} // namespace tm1p
