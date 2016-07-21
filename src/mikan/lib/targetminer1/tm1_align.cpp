#include <mikan/lib/targetminer1/include/tm1_align.hpp>         // TM1Alignment
#include <mikan/lib/targetminer1/include/tm1_inst_template.hpp> // TRNATYPE
#include <cmath>                 // max, min
#include <sstream>               // stringstream
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

namespace tm1p{

//
// TM1Alignment methods
//
template <class TRNAString>
void TM1Alignment<TRNAString>::clear_alignments()
{
    clear(mAlignBars);
    clear(mAlignMRNA);
    clear(mAlignMiRNA);
}

template <class TRNAString>
void TM1Alignment<TRNAString>::resize_alignments(unsigned pSize)
{
    resize(mAlignBars, pSize);
    resize(mAlignMRNA, pSize);
    resize(mAlignMiRNA, pSize);
}

template <class TRNAString>
void TM1Alignment<TRNAString>::align_all(
        TRNAString const &pMiRNASeq,
        TRNASet const &pMRNASeqs,
        TM1SeedSites<TRNAString> const &pSeedSites,
        String<int> const &pA1Pos)
{
    unsigned miLen;
    int a1pos;
    int sitepos;
    CharString matchBars;

    const String<unsigned>& mRNAPos = pSeedSites.get_mrna_pos();

    resize_alignments(length(pSeedSites.mEffectiveSites));

    miLen = length(pMiRNASeq);

   for (unsigned i = 0; i < length(mRNAPos); ++i)
   {
       if (!pSeedSites.mEffectiveSites[i])
       {
           continue;
       }

       resize(mAlignBars[i], miLen);
       resize(mAlignMRNA[i], miLen);
       resize(mAlignMiRNA[i], miLen);

       a1pos = pA1Pos[i];

       for (unsigned j = 0 ; j < miLen; ++j)
       {
           mAlignMiRNA[i][j] = pMiRNASeq[miLen - j - 1];
       }

       for (unsigned j = 0 ; j < miLen; ++j)
       {
           sitepos = a1pos - j;

           if (sitepos < 0 || sitepos > (int)length(pMRNASeqs[mRNAPos[i]]) - 1)
           {
               mAlignMRNA[i][miLen - j - 1] = ' ';
           }
           else
           {
               mAlignMRNA[i][miLen - j - 1] = pMRNASeqs[mRNAPos[i]][sitepos];
           }
       }

       for (unsigned j = 0 ; j < miLen; ++j)
       {
           if ((mAlignMiRNA[i][j] == 'A' && mAlignMRNA[i][j] == 'U')
                   || (mAlignMiRNA[i][j] == 'U' && mAlignMRNA[i][j] == 'A')
                   || (mAlignMiRNA[i][j] == 'G' && mAlignMRNA[i][j] == 'C')
                   || (mAlignMiRNA[i][j] == 'C' && mAlignMRNA[i][j] == 'G'))
           {
               mAlignBars[i][j] = '|';
           }
           else if ((mAlignMiRNA[i][j] == 'G' && mAlignMRNA[i][j] == 'U')
                   || (mAlignMiRNA[i][j] == 'U' && mAlignMRNA[i][j] == 'G'))
           {
               mAlignBars[i][j] = ':';
           }
           else
           {
               mAlignBars[i][j] = ' ';
           }

       }
   }

}

template <class TRNAString>
void TM1Alignment<TRNAString>::write_alignment(int pIdx) const
{
    std::stringstream stream;

    stream << "mRNA   5' " << mAlignMRNA[pIdx] << " 3'";
    stream << std::endl;
    stream << "          " << mAlignBars[pIdx] << "   ";
    stream << std::endl;
    stream << "miRNA  3' " << mAlignMiRNA[pIdx] << " 5'";
    stream << std::endl;

    std::cout << stream.str();
}

// Explicit template instantiation
template class TM1Alignment<TRNATYPE>;

} // namespace tm1p
