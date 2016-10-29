#include <tssvm_inst_template.hpp>  // TRNATYPE
#include <tssvm_seed_site.hpp>      // TSSVMSequences, TSSVMSeedSeqs, TSSVMSeedSites, TSSVMSeedSiteOverlap
#include <seqan/seq_io.h>

using namespace seqan;

namespace tssvm{

//
// TSSVMSequences methods
//
template <class TRNAString>
int TSSVMSequences<TRNAString>::read_fasta(CharString const &pFasta)
{
    CharString id;
    CharString seq;

    SequenceStream seqStream(toCString(pFasta));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file!" << std::endl;
        return 1;
    }

    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from " << toCString(pFasta) << "!" << std::endl;
            return 1;
        }

        toUpper(seq);
        for (unsigned i = 0; i < length(seq); ++i)
        {
            if (seq[i] == 'T')
            {
                seq[i] = 'U';
            }
        }
        appendValue(mSeqIds, id);
        appendValue(mSeqs, seq);

        if (mMaxLen < (int)length(seq))
        {
            mMaxLen = length(seq);
        }
    }

    return 0;
}

//
// TSSVMSeedSeqs methods
//
template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::create_seed_seqs()
{
    if (length(mMiRNASeq) == 0)
    {
        return 1;
    }

    TRNAString seedSeq;

    int retVal;

    resize(seedSeq, 6);
    for (unsigned i = 0; i < length(seedSeq); ++i)
    {
        seedSeq[i] = mMiRNASeq[i+1];
    }
    reverseComplement(seedSeq);

    appendValue(mSeedSeqs, seedSeq);
    appendValue(mSeedTypes, "6mer");
    appendValue(mMisMatchPos, 0);

    retVal = create_non_stringent_seed_seqs(seedSeq);
    if (retVal != 0)
    {
        return 1;
    }

    resize(mEffectiveSeeds, length(mSeedSeqs), true);

    return 0;
}

template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::create_non_stringent_seed_seqs(TRNAString& pSeedSeq)
{
    int retVal;

    retVal = create_guwobble_seed_seqs(pSeedSeq);
    if (retVal != 0)
    {
        return 1;
    }

    retVal = create_bt_seed_seqs(pSeedSeq);
    if (retVal != 0)
    {
        return 1;
    }

    retVal = create_bm_seed_seqs(pSeedSeq);
    if (retVal != 0)
    {
        return 1;
    }

    retVal = create_lp_seed_seqs(pSeedSeq);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::create_guwobble_seed_seqs(TRNAString& pSeedSeq)
{
    TRNAString seedGUSeq;
    TRNASet seedSeqs;
    StringSet<CharString> seedTypes;
    String<unsigned> misMatchPos;
    int retVal;

    // Add GUM first
    for (unsigned i = 0; i < length(pSeedSeq); ++i)
    {
        if (pSeedSeq[i] == 'C')
        {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            appendValue(seedSeqs, seedGUSeq);
            appendValue(seedTypes, "GUM");
            appendValue(misMatchPos, length(pSeedSeq)-i);
        }
        else if (pSeedSeq[i] == 'A')
        {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            appendValue(seedSeqs, seedGUSeq);
            appendValue(seedTypes, "GUT");
            appendValue(misMatchPos, length(pSeedSeq)-i);
        }

    }

    retVal = add_seeds_in_reverse_order(seedSeqs, seedTypes, misMatchPos);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::create_lp_seed_seqs(TRNAString& pSeedSeq)
{
    TRNAString seedLPSeq;
    bool effective;
    TRNASet seedSeqs;
    StringSet<CharString> seedTypes;
    String<unsigned> misMatchPos;
    int retVal;

    resize(seedLPSeq, 6);

    for (unsigned i = 0; i < length(pSeedSeq); ++i)
    {
        for (unsigned j = 0; j < length(mRNAChar); ++j)
        {
            effective = true;
            for (unsigned k = 0; k < length(pSeedSeq); ++k)
            {
                if (i == k)
                {
                    if ((mRNAChar[j] == pSeedSeq[k])
                            || (mRNAChar[j] ==  'G' && pSeedSeq[k] == 'A')
                            || (mRNAChar[j] ==  'U' && pSeedSeq[k] == 'C'))
                    {
                        effective = false;
                        break;
                    }

                    seedLPSeq[k] = mRNAChar[j];
                }
                else
                {
                    seedLPSeq[k] = pSeedSeq[k];
                }
            }

            if (effective)
            {
                appendValue(seedSeqs, seedLPSeq);
                appendValue(seedTypes, "LP");
                appendValue(misMatchPos, length(pSeedSeq)-i);
            }
        }
    }

    retVal = add_seeds_in_reverse_order(seedSeqs, seedTypes, misMatchPos);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::create_bm_seed_seqs(TRNAString& )
{
    TRNAString seedBMSeq;
    TRNASet seedSeqs;
    StringSet<CharString> seedTypes;
    String<unsigned> misMatchPos;
    int retVal;

    resize(seedBMSeq, 6);
    for (unsigned i = 0; i < length(seedBMSeq); ++i)
    {
        int k = 0;
        for (unsigned j = 0; j < length(seedBMSeq)+1; ++j)
        {
            if (i == j)
            {
                continue;
            }
            seedBMSeq[k] = mMiRNASeq[j+1];
            ++k;
        }
        reverseComplement(seedBMSeq);
        appendValue(seedSeqs, seedBMSeq);
        appendValue(seedTypes, "BM");
        appendValue(misMatchPos, i+1);
    }

    retVal = add_seeds_in_reverse_order(seedSeqs, seedTypes, misMatchPos);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::create_bt_seed_seqs(TRNAString& )
{
    TRNAString seedBTSeq;
    TRNASet seedSeqs;
    StringSet<CharString> seedTypes;
    String<unsigned> misMatchPos;
    int retVal;

    resize(seedBTSeq, 6);
    for (unsigned i = 0; i < length(seedBTSeq); ++i)
    {
        for (unsigned j = 0; j < length(mRNAChar); ++j)
        {
            int l = 0;
            for (unsigned k = 0; k < length(seedBTSeq); ++k)
            {
                if (i == k)
                {
                    seedBTSeq[k] = mRNAChar[j];
                }
                else
                {
                    seedBTSeq[k] = mMiRNASeq[l+2];
                    ++l;
                }
            }
            reverseComplement(seedBTSeq);
            appendValue(seedSeqs, seedBTSeq);
            appendValue(seedTypes, "BT");
            appendValue(misMatchPos, i+2);
        }
    }

    retVal = add_seeds_in_reverse_order(seedSeqs, seedTypes, misMatchPos);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

template <class TRNAString>
int TSSVMSeedSeqs<TRNAString>::add_seeds_in_reverse_order(
        TRNASet &pSeedSeqs,
        StringSet<CharString> &pSeedTypes,
        seqan::String<unsigned> &pMisMatchPos)
{
    resize(mSeedSeqs, length(mSeedSeqs)+length(pSeedSeqs));
    for (unsigned i = 0; i < length(pSeedSeqs); ++i)
    {
        mSeedSeqs[length(mSeedSeqs)-length(pSeedSeqs)+i] = pSeedSeqs[length(pSeedSeqs)-i-1];
    }

    resize(mSeedTypes, length(mSeedTypes)+length(pSeedTypes));
    for (unsigned i = 0; i < length(pSeedTypes); ++i)
    {
        mSeedTypes[length(mSeedTypes)-length(pSeedTypes)+i] = pSeedTypes[length(pSeedTypes)-i-1];
    }

    resize(mMisMatchPos, length(mMisMatchPos)+length(pMisMatchPos));
    for (unsigned i = 0; i < length(pMisMatchPos); ++i)
    {
        mMisMatchPos[length(mMisMatchPos)-length(pMisMatchPos)+i] = pMisMatchPos[length(pMisMatchPos)-i-1];
    }

    return 0;
}

template <class TRNAString>
void TSSVMSeedSeqs<TRNAString>::set_mirna_seq(TRNAString pSeq)
{
    clear(mSeedSeqs);
    clear(mMisMatchPos);
    clear(mEffectiveSeeds);
    mMiRNASeq = pSeq;
}

//
// TSSVMSeedSites methods
//
template <class TRNAString>
void TSSVMSeedSites<TRNAString>::reset_finder()
{
    goBegin(mFinder);
    clear(mFinder);
}

template <class TRNAString>
int TSSVMSeedSites<TRNAString>::find_seed_sites(
        TRNAString const &pMiRNA)
{
    TSSVMSeedSeqs<TRNAString> seedSeqs;
    TRNAString seedSeq;
    CharString newSeedType;
    int retVal;
    unsigned mRNAPos, sitePos;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0)
    {
        std::cerr << "ERROR: Could not get the seed sequence for " << pMiRNA;
        std::cerr << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(seedSeqs.mEffectiveSeeds); ++i)
    {
        if (!seedSeqs.mEffectiveSeeds[i])
        {
            continue;
        }

        seedSeq = seedSeqs.get_seed_seq(i);

        while (find(mFinder, seedSeq))
        {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            const CharString &seedType = seedSeqs.get_seed_type(i);
            if (seedType == "GUM" || seedType == "GUT" || seedType == "LP" || seedType == "BT")
            {
                if ((sitePos < MIN_DIST_TO_CDS+1) || (sitePos + MIN_DIST_UTR_END > length(mMRNASeqs[mRNAPos])-1))
                {
                    continue;
                }
            }
            else
            {
                if ((sitePos < MIN_DIST_TO_CDS) || (sitePos + MIN_DIST_UTR_END > length(mMRNASeqs[mRNAPos])))
                {
                    continue;
                }
            }

            retVal = set_seed_pos(seedSeqs, mRNAPos, sitePos, pMiRNA, seedSeq, i);
            if (retVal == 0)
            {
                appendValue(mMRNAPos, mRNAPos);
                appendValue(mSitePos, sitePos);
                appendValue(mEffectiveSites, true);
            }
        }
        reset_finder();
    }

    return 0;
}

template <class TRNAString>
void TSSVMSeedSites<TRNAString>::clear_pos()
{
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mS1Pos);
    clear(mS8Pos);
    clear(mMisMatchPos);
    clear(mEffectiveSites);
}

template <class TRNAString>
int TSSVMSeedSites<TRNAString>::set_seed_pos(
        TSSVMSeedSeqs<TRNAString> &pSeedSeqs,
        unsigned pMRNAPos,
        unsigned pSitePos,
        TRNAString const &pMiRNA,
        const TRNAString &pSeedSeq,
        unsigned pIdx)
{
    unsigned m8Pos, a1Pos;
    const CharString &seedType = pSeedSeqs.get_seed_type(pIdx);
    unsigned mmPos = pSeedSeqs.get_mis_matched_pos(pIdx);
    int retVal;

    m8Pos = pSitePos - 1;
    if (seedType == "BM")
    {
        m8Pos = pSitePos;
    }

    a1Pos = pSitePos + 6;
    if (seedType == "BT")
    {
        a1Pos = pSitePos + 7;
    }

    if ((seedType != "6mer") && (a1Pos >= length(mMRNASeqs[pMRNAPos])))
    {
        return 1;
    }

    retVal = set_seed_type(seedType, mMRNASeqs[pMRNAPos], pMiRNA, m8Pos, a1Pos, mmPos, pSeedSeq);
    if (retVal != 0)
    {
        return 1;
    }

    appendValue(mS1Pos, a1Pos);
    appendValue(mS8Pos, m8Pos);

    return 0;

}

template <class TRNAString>
int TSSVMSeedSites<TRNAString>::set_seed_type(
        const CharString &pCurType,
        const TRNAString &pMRNASeq,
        const TRNAString &pMiRNASeq,
        unsigned pM8Pos,
        unsigned pA1Pos,
        unsigned pMisMatchedPos,
        const TRNAString &)
{
    TRNAString newSeedSeq, complMiRNASeq;
    typename Value<TRNAString>::Type miRNAM2, miRNAM8;
    typename Value<TRNAString>::Type mRNAM2, mRNAM8, mRNAA1;
    CharString newSeedType;
    bool IsA1, MatchM8;

    complMiRNASeq = pMiRNASeq;
    complement(complMiRNASeq);
    miRNAM8 = complMiRNASeq[7];
    miRNAM2 = complMiRNASeq[1];

    mRNAM8 = pMRNASeq[pM8Pos];
    mRNAM2 = pMRNASeq[pA1Pos-1];
    if (pA1Pos < length(pMRNASeq))
    {
        mRNAA1 = pMRNASeq[pA1Pos];
    }
    else
    {
        mRNAA1 = 'A';
    }

    MatchM8 = false;
    if (miRNAM8 == mRNAM8)
    {
        MatchM8 = true;
    }

    IsA1 = false;
    if (mRNAA1 == 'A')
    {
        IsA1 = true;
    }

    newSeedType = "";
    if (pCurType == "6mer")
    {
        if (IsA1 && MatchM8)
        {
            newSeedType = "8mer";
        }
        else if (IsA1)
        {
            newSeedType = "7mer-A1";
        }
        else if (MatchM8)
        {
            newSeedType = "7mer-m8";
        }
        else
        {
            newSeedType = "6mer";
        }
    }
    else if (IsA1 && MatchM8)
    {
        if (pCurType == "BT")
        {
            if (miRNAM2 == mRNAM2)
            {
                newSeedType = pCurType;
            }
        }
        else
        {
            newSeedType = pCurType;
        }
    }

    if (newSeedType == "")
    {
        return 1;
    }

    appendValue(mSeedTypes, newSeedType);
    appendValue(mMisMatchPos, pMisMatchedPos);

    return 0;
}

//
// TSSVMSeedSiteOverlap methods
//

template <class TRNAString>
void TSSVMSeedSiteOverlap<TRNAString>::clear_site_pos()
{
    mRNAPosSet.clear();
    mSiteMap.clear();
    clear(mSortedMRNAPos);
}

template <class TRNAString>
void TSSVMSeedSiteOverlap<TRNAString>::cluster_site_pos(TSSVMSeedSites<TRNAString> &pSeedSites)
{
    const String<unsigned>& mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i)
    {
        mRNAPosSet.insert(mRNAPos[i]);
        mSiteMap.insert(std::pair<unsigned, unsigned>(mRNAPos[i], i));
    }
}

template <class TRNAString>
int TSSVMSeedSiteOverlap<TRNAString>::filter_overlapped_sites(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        unsigned pMRNANum)
{
    TItSet itSet;

    cluster_site_pos(pSeedSites);
    resize(mSortedMRNAPos, pMRNANum);

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet)
    {
        sort_by_seed_type(pSeedSites, *itSet);
    }

    return 0;
}

template <class TRNAString>
void TSSVMSeedSiteOverlap<TRNAString>::sort_by_seed_type(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        int pPosIdx)
{
    TItMap itMap;
    TItRetPair ret;
    std::multimap<unsigned, unsigned> sortedPos;
    std::multimap<unsigned, unsigned> sortedSeeds;
    StringSet<CharString> const& seedTypes = pSeedSites.get_seed_types();
    const String<unsigned>& sitePos = pSeedSites.get_site_pos();
    unsigned preced;
    unsigned count = 0;
    TITPos itPos;

    ret = mSiteMap.equal_range(pPosIdx);
    for (itMap = ret.first; itMap != ret.second; ++itMap)
    {
        sortedPos.insert(TPosPair((unsigned)sitePos[(*itMap).second], (*itMap).second));
        ++count;
    }

    resize(mSortedMRNAPos[pPosIdx], count);
    int n = 0;
    for (itPos = sortedPos.begin(); itPos != sortedPos.end(); ++itPos)
    {
        preced = get_seedtype_precedence(seedTypes[(*itPos).second]);
        sortedSeeds.insert(TPosPair(preced * (count + 1) + n, (*itPos).second));
        mSortedMRNAPos[pPosIdx][n] = (*itPos).second;
        ++n;
    }

    if (count > 1)
    {
        mark_overlapped_sites(pSeedSites, sortedSeeds);
    }
}

template <class TRNAString>
unsigned TSSVMSeedSiteOverlap<TRNAString>::get_seedtype_precedence(const CharString &pSeedType)
{
    unsigned preced;

    if (pSeedType == "8mer" || pSeedType == "7mer-A1" || pSeedType == "7mer-m8")
    {
        preced = 0;
    }
    else if (pSeedType == "6mer")
    {
        preced = 1;
    }
    else if (pSeedType == "GUM")
    {
        preced = 2;
    }
    else if (pSeedType == "GUT")
    {
        preced = 3;
    }
    else if (pSeedType == "BT")
    {
        preced = 4;
    }
    else if (pSeedType == "BM")
    {
        preced = 5;
    }
    else if (pSeedType == "LP")
    {
        preced = 6;
    }
    else
    {
        preced = 7;
    }

    return preced;
}

template <class TRNAString>
void TSSVMSeedSiteOverlap<TRNAString>::mark_overlapped_sites(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        std::multimap<unsigned, unsigned> &pSortedSeeds)
{
    IntervalTree<unsigned> tree;
    String<unsigned> results;
    const String<unsigned>& sitePos = pSeedSites.get_site_pos();
    StringSet<CharString> const& seedTypes = pSeedSites.get_seed_types();
    TITSeedTypes itSeedType;
    unsigned pos7, pos1, startPos, endPos;
    CharString seedType;
    bool searchOverlap;

    startPos = 0;
    endPos = 0;
    for (itSeedType = pSortedSeeds.begin(); itSeedType != pSortedSeeds.end(); ++itSeedType)
    {
        seedType = seedTypes[(*itSeedType).second];

        pos7 = sitePos[(*itSeedType).second] - 1;
        pos1 = sitePos[(*itSeedType).second] + 7;
        searchOverlap = true;

        if (seedType == "8mer" || seedType == "7mer-A1" || seedType == "7mer-m8" || seedType == "6mer")
        {
            startPos = pos7;
            endPos = pos1;
            searchOverlap = false;
        }
        else if (seedType == "GUM" || seedType == "GUT" || seedType == "LP")
        {
            startPos = pos7 - 1;
            endPos = pos1 + 1;
        }
        else if (seedType == "BT")
        {
            startPos = pos7 - 1;
            endPos = pos1 + 2;
        }
        else if (seedType == "BM")
        {
            pos7 = sitePos[(*itSeedType).second];
            pos1 = sitePos[(*itSeedType).second] + 8;
            startPos = pos7 - 1;
            endPos = pos1;
        }

        clear(results);
        if (searchOverlap)
        {
            findIntervals(tree, startPos, endPos, results);
        }

        if (length(results) == 0)
        {
            addInterval(tree, pos7, pos1);
        }

        if (searchOverlap && length(results) > 0)
        {
            pSeedSites.mEffectiveSites[(*itSeedType).second] = false;
        }
    }
}

// Explicit template instantiation
template class TSSVMSequences<TRNATYPE>;
template class TSSVMSeedSeqs<TRNATYPE>;
template class TSSVMSeedSites<TRNATYPE>;
template class TSSVMSeedSiteOverlap<TRNATYPE>;

} // namespace tssvm
