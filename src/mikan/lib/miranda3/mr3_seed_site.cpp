#include <mr3_inst_template.hpp>  // TRNATYPE
#include <mr3_seed_site.hpp>      // MR3Sequences, MR3SeedSeqs, MR3SeedSites
#include <seqan/seq_io.h>

using namespace seqan;

namespace mr3as{

//
// MR3Sequences methods
//
template <class TRNAString>
int MR3Sequences<TRNAString>::read_fasta(CharString const &pFasta)
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
            mMaxLen = (int)length(seq);
        }
    }

    return 0;
}

//
// MR3SeedSeqs methods
//
template <class TRNAString>
int MR3SeedSeqs<TRNAString>::create_seed_seqs(StringSet<CharString> &pSeedDef)
{
    if (length(mMiRNASeq) == 0)
    {
        return 1;
    }

    TRNAString seedSeq;
    CharString seedType;

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

    if (pSeedDef[2] == 'Y' || pSeedDef[1] == 'Y')
    {
        if (pSeedDef[3] == '1' || pSeedDef[3] == '+')
        {
            retVal = create_single_guwobble_seed_seqs(seedSeq);
        }
        if (pSeedDef[3] == '+')
        {
            retVal = create_multi_guwobble_seed_seqs(seedSeq);
        }

        if (pSeedDef[4] != "0:0")
        {
            retVal = create_mismatch_seed_seqs(seedSeq);
        }

        if (pSeedDef[4] != "0:0" && (pSeedDef[3] == '1' || pSeedDef[3] == '+'))
        {
            retVal = create_gu_mismatch_seed_seqs(seedSeq);
        }

        if (pSeedDef[5] == "1")
        {
            retVal = create_bt_seed_seqs(seedSeq);
        }

    }

    resize(mEffectiveSeeds, length(mSeedSeqs), true);
    retVal = check_redundant_seeds();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

template <class TRNAString>
int MR3SeedSeqs<TRNAString>::create_single_guwobble_seed_seqs(TRNAString &pSeedSeq)
{
    TRNAString seedGUSeq;

    for (unsigned i = 0; i < length(pSeedSeq); ++i)
    {
        if (pSeedSeq[i] == 'C')
        {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            appendValue(mSeedSeqs, seedGUSeq);
            appendValue(mSeedTypes, "GUT");
            appendValue(mMisMatchPos, i);
        }
        else if (pSeedSeq[i] == 'A')
        {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            appendValue(mSeedSeqs, seedGUSeq);
            appendValue(mSeedTypes, "GUM");
            appendValue(mMisMatchPos, i);
        }
    }

    return 0;
}

template <class TRNAString>
int MR3SeedSeqs<TRNAString>::create_multi_guwobble_seed_seqs(TRNAString &pSeedSeq)
{
    TRNAString seedGUSeq;

    unsigned seedDatLen = length(mSeedSeqs);

    for (unsigned i = 0; i < length(pSeedSeq); ++i)
    {
        if (pSeedSeq[i] == 'C')
        {
            for (unsigned j = 0; j < seedDatLen; ++j)
            {
                seedGUSeq = get_seed_seq(j);
                seedGUSeq[i] = 'U';
                appendValue(mSeedSeqs, seedGUSeq);
                appendValue(mSeedTypes, "GU+");
                appendValue(mMisMatchPos, 0);
            }
            seedDatLen = length(mSeedSeqs);
        }
        else if (pSeedSeq[i] == 'A')
        {
            for (unsigned j = 0; j < seedDatLen; ++j)
            {
                seedGUSeq = get_seed_seq(j);
                seedGUSeq[i] = 'G';
                appendValue(mSeedSeqs, seedGUSeq);
                appendValue(mSeedTypes, "GU+");
                appendValue(mMisMatchPos, 0);
            }
            seedDatLen = length(mSeedSeqs);
        }
    }

    return 0;
}

template <class TRNAString>
int MR3SeedSeqs<TRNAString>::create_mismatch_seed_seqs(TRNAString &pSeedSeq, bool pIsMMGU, int pGUPos)
{
    TRNAString seedLPSeq;
    CharString seedType;
    char ch1=0;
    char ch2=0;
    char ch3=0;

    if (pIsMMGU)
    {
        seedType ="MMGU";
    }
    else
    {
        seedType ="MM";
    }

    for (unsigned i = 0; i < length(pSeedSeq); ++i)
    {
        if (pIsMMGU && (pGUPos == (int)i))
        {
            continue;
        }
        if (pSeedSeq[i] == 'A')
        {
            ch1 = 'C';
            ch2 = 'x';
            ch3 = 'U';
        }
        else if (pSeedSeq[i] == 'C')
        {
            ch1 = 'A';
            ch2 = 'x';
            ch3 = 'G';
        }
        else if (pSeedSeq[i] == 'G')
        {
            ch1 = 'A';
            ch2 = 'U';
            ch3 = 'C';
        }
        else if (pSeedSeq[i] == 'U')
        {
            ch1 = 'C';
            ch2 = 'G';
            ch3 = 'A';
        }

        seedLPSeq = pSeedSeq;
        seedLPSeq[i] = ch1;
        appendValue(mSeedSeqs, seedLPSeq);
        appendValue(mSeedTypes, seedType);
        appendValue(mMisMatchPos, i);
        if (ch2 != 'x')
        {
            seedLPSeq = pSeedSeq;
            seedLPSeq[i] = ch2;
            appendValue(mSeedSeqs, seedLPSeq);
            appendValue(mSeedTypes, seedType);
            appendValue(mMisMatchPos, i);
        }
        seedLPSeq = pSeedSeq;
        seedLPSeq[i] = ch3;
        appendValue(mSeedSeqs, seedLPSeq);
        appendValue(mSeedTypes, seedType);
        appendValue(mMisMatchPos, i);
    }

    return 0;
}

template <class TRNAString>
int MR3SeedSeqs<TRNAString>::create_gu_mismatch_seed_seqs(TRNAString &pSeedSeq)
{
    TRNAString seedGUSeq;

    for (unsigned i = 0; i < length(pSeedSeq); ++i)
    {
        if (pSeedSeq[i] == 'C')
        {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            create_mismatch_seed_seqs(seedGUSeq, true, i);
        }
        else if (pSeedSeq[i] == 'A')
        {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            create_mismatch_seed_seqs(seedGUSeq, true, i);
        }
    }

    return 0;
}

template <class TRNAString>
int MR3SeedSeqs<TRNAString>::create_bt_seed_seqs(TRNAString &pSeedSeq)
{
    TRNAString seedBTSeq;
    TRNAString allRNAs;

    resize(allRNAs, 4);
    allRNAs[0] = 'A';
    allRNAs[1] = 'C';
    allRNAs[2] = 'G';
    allRNAs[3] = 'U';

    resize(seedBTSeq, 6);
    for (unsigned i = 0; i < length(seedBTSeq); ++i)
    {
        for (unsigned j = 0; j < length(allRNAs); ++j)
        {
            int l = 0;
            for (unsigned k = 0; k < length(seedBTSeq); ++k)
            {
                if (i == k)
                {
                    seedBTSeq[k] = allRNAs[j];
                    appendValue(mMisMatchPos, i);
                }
                else
                {
                    seedBTSeq[k] = pSeedSeq[l];
                    ++l;
                }
            }
            appendValue(mSeedSeqs, seedBTSeq);
            appendValue(mSeedTypes, "BT");
        }
    }

    return 0;
}


template <class TRNAString>
int MR3SeedSeqs<TRNAString>::check_redundant_seeds()
{
    typedef Index<StringSet<TRNAString>, IndexQGram<UngappedShape<6> > > TIndexQGram;
    typedef Finder<TIndexQGram> TFinder;

    TRNAString seedSeq;
    TIndexQGram RNAIdx(mSeedSeqs);
    TFinder finder(RNAIdx);

    for (unsigned i = 0; i < length(mSeedSeqs); ++i)
    {
        if (!mEffectiveSeeds[i])
        {
            continue;
        }
        seedSeq = mSeedSeqs[i];
        while (find(finder, seedSeq))
        {
            if (i != position(finder).i1)
            {
                mEffectiveSeeds[position(finder).i1] = false;
            }
        }
        goBegin(finder);
        clear(finder);
    }

    return 0;
}

template <class TRNAString>
void MR3SeedSeqs<TRNAString>::set_mirna_seq(TRNAString pSeq)
{
    clear(mSeedSeqs);
    clear(mSeedTypes);
    clear(mEffectiveSeeds);
    mMiRNASeq = pSeq;
}

//
// MR3SeedSites methods
//
template <class TRNAString>
void MR3SeedSites<TRNAString>::clear_pos()
{
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mMisMatchPos);
    clear(mEffectiveSites);
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::reset_finder()
{
    goBegin(mFinder);
    clear(mFinder);
}

template <class TRNAString>
int MR3SeedSites<TRNAString>::find_seed_sites(
        TRNAString const &pMiRNA,
        StringSet<CharString> &pSeedDef)
{
    MR3SeedSeqs<TRNAString> seedSeqs;
    TRNAString seedSeq;
    CharString seedType;
    int retVal;
    unsigned mRNAPos, sitePos;
    bool effectiveSite;
    unsigned misMatchPos;
    unsigned endPos;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    retVal = seedSeqs.create_seed_seqs(pSeedDef);
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
        seedType = seedSeqs.get_seed_type(i);
        misMatchPos = seedSeqs.get_mismatched_pos(i);

        while (find(mFinder, seedSeq))
        {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            effectiveSite = true;
            endPos = sitePos + INDEXED_SEQ_LEN;
            if ((endPos < MIN_DIST_TO_CDS) || (endPos + MIN_DIST_UTR_END > length(mMRNASeqs[mRNAPos])))
            {
                effectiveSite = false;
            }

            appendValue(mMRNAPos, mRNAPos);
            appendValue(mSitePos, sitePos);
            set_new_seed_type(seedType, pSeedDef, mRNAPos, sitePos, pMiRNA, misMatchPos, effectiveSite);
            appendValue(mEffectiveSites, effectiveSite);
        }
        reset_finder();
    }

    return 0;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_new_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        unsigned pMRNAPos,
        unsigned pSitePos,
        TRNAString const &pMiRNA,
        unsigned pMisMatchPos,
        bool &pEffectiveSite)
{
    bool matchM8, matchM9, gutM8, gutM9, gumM8, gumM9;;
    CharString newSeedType = "";

    if (!pEffectiveSite)
    {
        appendValue(mSeedTypes, "");
        appendValue(mMisMatchPos, 0);
        return;
    }

    set_mx_matches(pMRNAPos, pSitePos, pMiRNA, 8, matchM8, gutM8, gumM8);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNA, 9, matchM9, gutM9, gumM9);

    set_stringent_seed_type(pCurSeedType, pSeedDef, matchM8, matchM9, pMisMatchPos, newSeedType);

    if (newSeedType == "" && (pSeedDef[2] == 'Y' || pSeedDef[1] == 'Y'))
    {
        if (pSeedDef[3] == '1' || pSeedDef[3] == '+')
        {
            set_single_gu_seed_type(pCurSeedType, pSeedDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                    pMisMatchPos, newSeedType);
        }

        if (newSeedType == "" && pSeedDef[3] == '+')
        {
            set_multiple_gu_seed_type(pCurSeedType, pSeedDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                    pMisMatchPos, newSeedType);
        }

        if (newSeedType == "" && pSeedDef[4] != "0:0")
        {
            set_mismatch_seed_type(pCurSeedType, pSeedDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                    pMisMatchPos, newSeedType);
        }

        if (newSeedType == "" && pSeedDef[4] != "0:0" && (pSeedDef[3] == '1' || pSeedDef[3] == '+'))
        {
            set_gu_mismatch_seed_type(pCurSeedType, pSeedDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                    pMisMatchPos, newSeedType);
        }

        if (pCurSeedType == "BT" && newSeedType == "" && pSeedDef[5] == "1")
        {
            set_bt_seed_type(pMRNAPos, pSitePos, pMiRNA, pMisMatchPos, newSeedType);
        }

    }

    set_6mer_seed_type(pCurSeedType, pSeedDef, matchM8, matchM9, pMisMatchPos, newSeedType);

    if (newSeedType == "")
    {
        pEffectiveSite = false;
        appendValue(mSeedTypes, pCurSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;

}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_mx_matches(
        unsigned pMRNAPos,
        unsigned pSitePos,
        TRNAString const &pMiRNA,
        int pMx,
        bool &pMatchMx,
        bool &pGutMx,
        bool &pGumMx)
{
    TRNAString cMiRNASeq, miRNAMx, mRNAMx, miRNAMxC;

    miRNAMx = pMiRNA[pMx-1];
    cMiRNASeq = pMiRNA;
    complement(cMiRNASeq);
    miRNAMxC = cMiRNASeq[pMx-1];

    pMatchMx = false;
    pGutMx = false;
    pGumMx = false;

    if ((int)pSitePos - pMx + 1 + (int)INDEXED_SEQ_LEN < 0)
    {
        return;
    }
    mRNAMx = mMRNASeqs[pMRNAPos][pSitePos - (pMx - 1 - INDEXED_SEQ_LEN)];

    pMatchMx = false;
    if (miRNAMxC == mRNAMx)
    {
        pMatchMx = true;
    }

    pGutMx = false;
    pGumMx = false;
    if (miRNAMx == 'G' && mRNAMx == 'U')
    {
        pGutMx = true;
    }
    else if (miRNAMx == 'U' && mRNAMx == 'G')
    {
        pGumMx = true;
    }

}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_stringent_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        bool pMatchMx1,
        bool pMatchMx2,
        unsigned,
        CharString &pNewSeedType)
{
    if (pNewSeedType != "")
    {
        return;
    }

    if (pCurSeedType == "6mer")
    {
        if (pSeedDef[2] == 'Y' && pMatchMx1 && pMatchMx2)
        {
            pNewSeedType = "8mer";
        }
        else if (pSeedDef[1] == 'Y' && pMatchMx1)
        {
            pNewSeedType = "7mer";
        }
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_6mer_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        bool ,
        bool ,
        unsigned,
        CharString &pNewSeedType)
{
    if (pNewSeedType != "")
    {
        return;
    }

    if (pCurSeedType == "6mer")
    {
        if (pSeedDef[0] == 'Y')
        {
            pNewSeedType = "6mer";
        }
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_single_gu_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned pMisMatchPos,
        CharString &pNewSeedType)
{
    int mm;

    if (pNewSeedType != "")
    {
        return;
    }

    if (pCurSeedType == "6mer")
    {
        if (pSeedDef[2] == 'Y' && ((pMatchMx1 && pGutMx2) || (pMatchMx2 && pGutMx1)))
        {
            pNewSeedType = "8mer_GUT";
            if (pMatchMx1)
            {
                mm = pM2;
            }
            else
            {
                mm = pM1;
            }
        }
        else if (pSeedDef[2] == 'Y' && ((pMatchMx1 && pGumMx2) || (pMatchMx2 && pGumMx1)))
        {
            pNewSeedType = "8mer_GUM";
            if (pMatchMx1)
            {
                mm = pM2;
            }
            else
            {
                mm = pM1;
            }
        }
        else if (pSeedDef[1] == 'Y' && pGutMx1)
        {
            pNewSeedType = "7mer_GUT";
            mm = pM1;
        }
        else if (pSeedDef[1] == 'Y' && pGumMx1)
        {
            pNewSeedType = "7mer_GUM";
            mm = pM1;
        }
    }
    else if (pCurSeedType == "GUT")
    {
        if (pSeedDef[2] == 'Y' && pMatchMx1 && pMatchMx2)
        {
            pNewSeedType = "8mer_GUT";
            mm = pMisMatchPos;
        }
        else if (pSeedDef[1] == 'Y' && pMatchMx1)
        {
            pNewSeedType = "7mer_GUT";
            mm = pMisMatchPos;
        }
    }
    else if (pCurSeedType == "GUM")
    {
        if (pSeedDef[2] == 'Y' && pMatchMx1 && pMatchMx2)
        {
            pNewSeedType = "8mer_GUM";
            mm = pMisMatchPos;
        }
        else if (pSeedDef[1] == 'Y' && pMatchMx1)
        {
            pNewSeedType = "7mer_GUM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "")
    {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_multiple_gu_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int,
        int,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned,
        CharString &pNewSeedType)
{
    if (pNewSeedType != "")
    {
        return;
    }

    if (pCurSeedType == "6mer" && pSeedDef[2] == 'Y' && (pGumMx1 || pGutMx1) && (pGumMx2 || pGutMx2))
    {
        pNewSeedType = "8mer_GU+";
    }
    else if (pCurSeedType == "GUT" || pCurSeedType == "GUM")
    {
        if (pSeedDef[2] == 'Y' && (pGumMx1 || pGutMx1 || pMatchMx1) && (pGumMx2 || pGutMx2 || pMatchMx2))
        {
            pNewSeedType = "8mer_GU+";
        }
        else if (pSeedDef[1] == 'Y' && (pGumMx1 || pGutMx1))
        {
            pNewSeedType = "7mer_GU+";
        }
    }
    else if (pCurSeedType == "GU+")
    {
        if (pSeedDef[2] == 'Y' && (pGumMx1 || pGutMx1 || pMatchMx1) && (pGumMx2 || pGutMx2 || pMatchMx2))
        {
            pNewSeedType = "8mer_GU+";
        }
        else if (pSeedDef[1] == 'Y' && (pGumMx1 || pGutMx1 || pMatchMx1))
        {
            pNewSeedType = "7mer_GU+";
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "")
    {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_mismatch_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned pMisMatchPos,
        CharString &pNewSeedType)
{
    int mm;

    if (pNewSeedType != "")
    {
        return;
    }

    if (pSeedDef[2] == 'Y')
    {
        if (pCurSeedType == "6mer" && ((pMatchMx1 && !pMatchMx2 && !pGutMx2 && !pGumMx2)
                || (!pMatchMx1 && !pGutMx1 && !pGumMx1 && pMatchMx2)))
        {
            pNewSeedType = "8mer_MM";
            if (pMatchMx1)
            {
                mm = pM2;
            }
            else
            {
                mm = pM1;
            }
        }
        else if (pCurSeedType == "MM" && pMatchMx1 && pMatchMx2)
        {
            pNewSeedType = "8mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1")
    {
        if (pCurSeedType == "6mer" && !pMatchMx1 && !pGutMx1 && !pGumMx1)
        {
            pNewSeedType = "7mer_MM";
            mm = pM1;
        }
        else if (pCurSeedType == "MM" && pMatchMx1)
        {
            pNewSeedType = "7mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "")
    {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_gu_mismatch_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned pMisMatchPos,
        CharString &pNewSeedType)
{
    int mm;

    if (pNewSeedType != "")
    {
        return;
    }

    if (pSeedDef[2] == 'Y')
    {
        if (pCurSeedType == "6mer" && ((!pMatchMx1 && !pMatchMx2 && (pGutMx1 || pGumMx1) && !(pGutMx2 || pGumMx2))
                || (!pMatchMx1 && !pMatchMx2 && !(pGutMx1 || pGumMx1) && (pGutMx2 || pGumMx2))))
        {
            pNewSeedType = "8mer_MMGU";
            if (pGutMx1 || pGumMx1)
            {
                mm = pM2;
            }
            else
            {
                mm = pM1;
            }
        }
        else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") &&
                ((pMatchMx1 && !pMatchMx2 && !pGutMx2 && !pGumMx2)
                || (!pMatchMx1 && pMatchMx2 && !pGutMx1 && !pGumMx1)))
        {
            pNewSeedType = "8mer_MMGU";
            if (pMatchMx1)
            {
                mm = pM2;
            }
            else
            {
                mm = pM1;
            }
        }
        else if (pCurSeedType == "MM" && ((pMatchMx1 && (pGutMx2 || pGumMx2))
                || ((pGutMx1 || pGumMx1) && pMatchMx2)))
        {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        }
        else if (pCurSeedType == "MMGU" && pMatchMx1 && pMatchMx2)
        {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1")
    {
        if (pCurSeedType == "MM" && (pGutMx2 || pGumMx2))
        {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        }
        else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") && (!pMatchMx1 && !pGutMx1 && !pGumMx1))
        {
            pNewSeedType = "7mer_MMGU";
            mm = pM1;
        }
        else if (pCurSeedType == "MMGU" && pMatchMx1)
        {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "")
    {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

template <class TRNAString>
void MR3SeedSites<TRNAString>::set_bt_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        TRNAString const &pMiRNA,
        unsigned pMisMatchPos,
        CharString &pNewSeedType)
{
    TRNAString cMiRNASeq, miRNAM6, miRNAM6C, miRNAM7, miRNAM7C, miRNAM8, miRNAM8C, mRNAM7, mRNAM8, mRNAM9;

    if (pNewSeedType != "")
    {
        return;
    }

    cMiRNASeq = pMiRNA;
    complement(cMiRNASeq);
    miRNAM6 = pMiRNA[6];
    miRNAM6C = cMiRNASeq[6];
    miRNAM7 = pMiRNA[7];
    miRNAM7C = cMiRNASeq[7];
    miRNAM8 = pMiRNA[8];
    miRNAM8C = cMiRNASeq[8];

    if ((int)pSitePos - 7 + (int)INDEXED_SEQ_LEN < 0)
    {
        return;
    }
    mRNAM7 = mMRNASeqs[pMRNAPos][pSitePos- (7 - INDEXED_SEQ_LEN)];

    if ((int)pSitePos - 8 + (int)INDEXED_SEQ_LEN < 0)
    {
        return;
    }
    mRNAM8 = mMRNASeqs[pMRNAPos][pSitePos- (8 - INDEXED_SEQ_LEN)];

    if (miRNAM6C == mRNAM7 && miRNAM7C == mRNAM8)
    {
        pNewSeedType = "7mer_BT";
    }

    if (pNewSeedType == "7mer_BT" && ((int)pSitePos - 9 + (int)INDEXED_SEQ_LEN >= 0))
    {
        mRNAM9 = mMRNASeqs[pMRNAPos][pSitePos- (9 - INDEXED_SEQ_LEN)];
        if (miRNAM8C == mRNAM9)
        {
            pNewSeedType = "8mer_BT";
        }
    }

    //TODO Check BT sites - seed pos 3-7 instead of pos 2-6?
//    std::cout << std::endl;
//    std::cout << mMRNASeqs[pMRNAPos][pSitePos] << "," <<  mMRNASeqs[pMRNAPos][pSitePos+1] << "," <<  mMRNASeqs[pMRNAPos][pSitePos+2] << std::endl;
//    std::cout << mRNAM7 << "," <<  mRNAM8 << std::endl;
//    std::cout << miRNAM6C << "," <<  miRNAM7C << std::endl;

    if (FORCE_LAST_MATCH && pNewSeedType != "8mer_BT")
    {
        pNewSeedType = "";
    }

    if (pNewSeedType != "")
    {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, pMisMatchPos);
    }

    return;
}


template <class TRNAString>
void MR3SeedSites<TRNAString>::check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType)
{
    if (pNewSeedType == "6mer" || pNewSeedType == "7mer" || pNewSeedType == "8mer")
    {
        return;
    }

    if ((pNewSeedType[0] == '7' && !pMatchM8) || (pNewSeedType[0] == '8' && !pMatchM9))
    {
        pNewSeedType = "";
    }
}

// Explicit template instantiation
template class MR3Sequences<TRNATYPE>;
template class MR3SeedSeqs<TRNATYPE>;
template class MR3SeedSites<TRNATYPE>;

} // namespace mr3as
