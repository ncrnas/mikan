#include <tssvm_align.hpp>          // TSAlign
#include <tssvm_inst_template.hpp>  // TRNATYPE

using namespace seqan;

namespace tssvm{

//
// TSAlign methods
//
template <class TRNAString>
void TSAlign<TRNAString>::clear_alignments()
{
    clear(mAlignBars);
    clear(mAlignMRNA);
    clear(mAlignMiRNA);
    clear(mAlignScores);
}

template <class TRNAString>
int TSAlign<TRNAString>::align_seq(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TRNAString const &pMiRNASeq,
        TRNASet const &pMRNASeqs)
{
    const String<unsigned>& mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned>& sitePos = pSeedSites.get_site_pos();
    const StringSet<CharString>& seedTypes = pSeedSites.get_seed_types();
    const String<unsigned>& misMatchPos = pSeedSites.get_mismatch_pos();
    TRNAString miRNAAlignSeq;
    TRNAString mRNAAlignSeq;
    TAlign align;
    unsigned alignLen;
    int retVal;

    resize(mAlignBars, length(pSeedSites.mEffectiveSites));
    resize(mAlignMRNA, length(pSeedSites.mEffectiveSites));
    resize(mAlignMiRNA, length(pSeedSites.mEffectiveSites));
    resize(mAlignScores, length(pSeedSites.mEffectiveSites));

    for (unsigned i = 0; i < length(mRNAPos); ++i)
    {
        if (!pSeedSites.mEffectiveSites[i])
        {
            continue;
        }

        retVal = set_mirna_seq_for_align(seedTypes[i], pMiRNASeq, miRNAAlignSeq);
        if (retVal != 0)
        {
            return 1;
        }

        retVal = set_mrna_seq_for_align(seedTypes[i], (unsigned)sitePos[i], length(pMiRNASeq),
                pMRNASeqs[mRNAPos[i]], mRNAAlignSeq);
        if (retVal != 0)
        {
            return 1;
        }

        retVal = set_addtional_sequences(miRNAAlignSeq, mRNAAlignSeq);
        if (retVal != 0)
        {
            return 1;
        }

        clearClipping(align);
        clearGaps(align);
        resize(rows(align), 2);
        assignSource(row(align, 0), miRNAAlignSeq);
        assignSource(row(align, 1), mRNAAlignSeq);

        mAlignScores[i] = globalAlignment(align, mScoreMatrix, AlignConfig<true, true, false, false>());

        alignLen = get_align_len(align);

        retVal = set_align_mrna(align, alignLen, seedTypes[i], (unsigned)misMatchPos[i], (unsigned)sitePos[i],
                pMRNASeqs[mRNAPos[i]], i);
        if (retVal != 0)
        {
            return 1;
        }

        retVal = set_align_mirna(align, alignLen, seedTypes[i], (unsigned)misMatchPos[i], pMiRNASeq, i);
        if (retVal != 0)
        {
            return 1;
        }

        retVal = set_align_bars(i, alignLen);
        if (retVal != 0)
        {
            return 1;
        }

    }

    return 0;
}

template <class TRNAString>
int TSAlign<TRNAString>::set_mirna_seq_for_align(
        const CharString &pSeedType,
        TRNAString const &pMiRNASeq,
        TRNAString &pMiRNAAlignSeq)
{
    int seqStart, seqEnd;

    seqStart = 8;
    seqEnd = length(pMiRNASeq);
    if (pSeedType == "BT")
    {
        seqStart = 7;
    }

    resize(pMiRNAAlignSeq, seqEnd - seqStart);;
    for (int i = 0; i < seqEnd - seqStart; ++i)
    {
        pMiRNAAlignSeq[i] = pMiRNASeq[seqStart+i];
    }

    reverse(pMiRNAAlignSeq);

    return 0;

}

template <class TRNAString>
int TSAlign<TRNAString>::set_mrna_seq_for_align(
        const CharString &pSeedType,
        unsigned pSitePos,
        unsigned pMiRNALen,
        const TRNAString &pMRNASeq,
        TRNAString &pMRNAAlignSeq)
{
    int seqStart, seqEnd;

    if (pSeedType  == "BM")
    {
        seqEnd = pSitePos;
        seqStart = seqEnd - (pMiRNALen - 8);
    }
    else if (pSeedType  == "BT")
    {
        seqEnd = pSitePos;
        seqStart = seqEnd - (pMiRNALen - 7);
    }
    else
    {
        seqEnd = pSitePos - 1;
        seqStart = seqEnd - (pMiRNALen - 8);
    }

    resize(pMRNAAlignSeq, seqEnd - seqStart);
    if (seqStart < 0)
    {
        for (int i = 0; i < seqStart * -1; ++i)
        {
            appendValue(pMRNAAlignSeq, 'A');
        }
        seqStart = 0;
    }

    for (int i = 0; i < seqEnd - seqStart; ++i)
    {
        pMRNAAlignSeq[i] = pMRNASeq[seqStart+i];
    }

    return 0;

}

template <class TRNAString>
int TSAlign<TRNAString>::set_addtional_sequences(TRNAString &pMiRNAAlignSeq, TRNAString &pMRNAAlignSeq)
{
    typename Value<TRNAString>::Type miRNAEndChar;
    typename Value<TRNAString>::Type miAddChar;
    typename Value<TRNAString>::Type mRNAEndChar;
    typename Value<TRNAString>::Type mAddChar;

    miRNAEndChar = pMiRNAAlignSeq[length(pMiRNAAlignSeq)-1];
    mRNAEndChar = pMRNAAlignSeq[length(pMRNAAlignSeq)-1];

    if ((mRNAEndChar == 'A' || mRNAEndChar == 'U') && (miRNAEndChar == 'A' || miRNAEndChar == 'U'))
    {
        if (mRNAEndChar == 'A')
        {
            mAddChar = 'G';
            miAddChar = 'C';
        }
        else
        {
            mAddChar = 'C';
            miAddChar = 'G';
        }
    }
    else if ((mRNAEndChar == 'C' || mRNAEndChar == 'G') && (miRNAEndChar == 'C' || miRNAEndChar == 'G'))
    {
        if (mRNAEndChar == 'G')
        {
            mAddChar = 'U';
            miAddChar = 'A';
        }
        else
        {
            mAddChar = 'A';
            miAddChar = 'U';
        }
    }
    else
    {
        if (mRNAEndChar == 'C' || mRNAEndChar == 'G')
        {
            mAddChar = 'G';
            miAddChar = 'C';
        }
        else
        {
            mAddChar = 'C';
            miAddChar = 'G';
        }
    }

    resize(pMRNAAlignSeq, length(pMRNAAlignSeq)+8, mAddChar);
    resize(pMiRNAAlignSeq, length(pMiRNAAlignSeq)+8, miAddChar);

    resize(pMRNAAlignSeq, length(pMRNAAlignSeq)+8);
    resize(pMiRNAAlignSeq, length(pMiRNAAlignSeq)+8);

    replace(pMRNAAlignSeq, length(pMRNAAlignSeq)-8, length(pMRNAAlignSeq), "CGCGAUAU");
    replace(pMiRNAAlignSeq, length(pMiRNAAlignSeq)-8, length(pMiRNAAlignSeq), "GCGCUAUA");

    return 0;
}

template <class TRNAString>
unsigned TSAlign<TRNAString>::get_align_len(TAlign &pAlign)
{
    unsigned alignLen;
    TGap &miRNAAlign = row(pAlign, 0);
    TGap &mRNAAlign = row(pAlign, 1);

    if (length(miRNAAlign) > length(mRNAAlign))
    {
        alignLen = length(miRNAAlign);
    }
    else
    {
        alignLen = length(mRNAAlign);
    }

    alignLen -= 16;

    return alignLen;

}

template <class TRNAString>
int TSAlign<TRNAString>::set_align_mrna(
        TAlign &pAlign,
        unsigned pAlignLen,
        const CharString &pSeedType,
        unsigned pMisMatchPos,
        unsigned pSitePos,
        const TRNAString &pMRNASeq,
        int pIdx)
{
    TGap &mRNAAlign = row(pAlign, 1);
    CharString seedSiteSeq;
    int seqStart;
    unsigned seedSiteSeqLen;

    resize(mAlignMRNA[pIdx], pAlignLen+8);
    for (unsigned i = 0; i < pAlignLen; ++i)
    {
        mAlignMRNA[pIdx][i] = mRNAAlign[i];
    }

    seqStart = pSitePos - 1;
    if (pSeedType  == "BM")
    {
        seqStart = pSitePos;
    }

    seedSiteSeqLen = 8;
    if (pSeedType  == "BT")
    {
        seedSiteSeqLen = 9;
    }

    resize(seedSiteSeq, seedSiteSeqLen);
    unsigned k = 0;
    for (unsigned i = 0; i < length(seedSiteSeq); ++i)
    {
        if (seqStart+k < length(pMRNASeq))
        {
            if ((pSeedType == "BM") && (pMisMatchPos == length(seedSiteSeq) - (i+1)))
            {
                seedSiteSeq[i] = '-';
            }
            else
            {
                seedSiteSeq[i] = pMRNASeq[seqStart+k];
                ++k;
            }
        }
        else
        {
            seedSiteSeq[i] = 'A';
        }
    }
    replace(mAlignMRNA[pIdx], length(mAlignMRNA[pIdx])-seedSiteSeqLen, length(mAlignMRNA[pIdx]),
            toCString(seedSiteSeq));

    return 0;
}

template <class TRNAString>
int TSAlign<TRNAString>::set_align_mirna(
        TAlign &pAlign,
        unsigned pAlignLen,
        const CharString &pSeedType,
        unsigned pMisMatchPos,
        const TRNAString &pMiRNASeq,
        int pIdx)
{
    TGap &miRNAAlign = row(pAlign, 0);
    CharString seedSeq;

    resize(mAlignMiRNA[pIdx], pAlignLen+8);
    for (unsigned i = 0; i < pAlignLen; ++i)
    {
        mAlignMiRNA[pIdx][i] = miRNAAlign[i];
    }

    resize(seedSeq, 8);
    unsigned k = 0;
    for (unsigned i = 0; i < length(seedSeq); ++i)
    {
        if ((pSeedType == "BT") && (pMisMatchPos == i))
        {
            seedSeq[i] = '-';
        }
        else
        {
            seedSeq[i] = pMiRNASeq[k];
            ++k;
        }
    }
    reverse(seedSeq);
    replace(mAlignMiRNA[pIdx], length(mAlignMiRNA[pIdx])-8, length(mAlignMiRNA[pIdx]), toCString(seedSeq));

    return 0;
}

template <class TRNAString>
int TSAlign<TRNAString>::set_align_bars(int pIdx, unsigned pAlignLen)
{
    char miRNA, mRNA;

    resize(mAlignBars[pIdx], pAlignLen+8);

    for (unsigned i = 0; i < length(mAlignBars[pIdx]); ++i)
    {
        miRNA = mAlignMiRNA[pIdx][i];
        mRNA = mAlignMRNA[pIdx][i];

        if (miRNA ==  'A')
        {
            if (mRNA == 'U')
            {
                mAlignBars[pIdx][i] = '|';
            }
            else
            {
                mAlignBars[pIdx][i] = ' ';
            }
        }
        else if (miRNA ==  'C')
        {
            if (mRNA == 'G')
            {
                mAlignBars[pIdx][i] = '|';
            }
            else
            {
                mAlignBars[pIdx][i] = ' ';
            }
        }
        else if (miRNA ==  'G')
        {
            if (mRNA == 'C')
            {
                mAlignBars[pIdx][i] = '|';
            }
            else if (mRNA == 'U')
            {
                mAlignBars[pIdx][i] = ':';
            }
            else
            {
                mAlignBars[pIdx][i] = ' ';
            }
        }
        else if (miRNA ==  'U')
        {
            if (mRNA == 'A')
            {
                mAlignBars[pIdx][i] = '|';
            }
            else if (mRNA == 'G')
            {
                mAlignBars[pIdx][i] = ':';
            }
            else
            {
                mAlignBars[pIdx][i] = ' ';
            }
        }
        else
        {
            mAlignBars[pIdx][i] = ' ';
        }
    }

    return 0;
}

template <class TRNAString>
void TSAlign<TRNAString>::write_alignment(int pIdx)
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
template class TSAlign<TRNATYPE>;

} // namespace tssvm
