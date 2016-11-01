#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif
#include <seqan/arg_parse.h>
#include <mr3_inst_template.hpp>  // TRNATYPE
#include <mr3_option.hpp>         // MR3Options
#include <mr3_seed_site.hpp>      // MR3Sequences, MR3SeedSites
#include <mr3_score.hpp>          // MR3SiteScores, MR3TotalScores
#include <mr3_site_cluster.hpp>   // MR3Overlap, MR3TopNScore, MR3SortedSitePos
#include <mr3_core.hpp>           // MR3CoreInput, MR3Core


namespace mr3as {

int MR3CoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line.
    mr3as::MR3Options options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return retVal;
    }

    // Create index
    mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    mr3as::MR3Core<mr3as::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();

    mr3as::MR3Core<mr3as::TRNATYPE> pitaCore(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    pitaCore.init_from_args(options);
    pitaCore.open_output_file();
    retVal = pitaCore.calculate_all_scores();

    return retVal;
}

//
// MR3CoreInput methods
//
template <class TRNAString>
void MR3CoreInput<TRNAString>::init_from_args(MR3Options& opts)
{
    mMiRNAFasta = opts.mMiRNAFasta;
    mMRNAFasta = opts.mMRNAFasta;
}

template <class TRNAString>
int MR3CoreInput<TRNAString>::load_seq_from_file()
{
    int retVal;

    // Read miRNA fasta file
    retVal = mMiRNASeqs.read_fasta(mMiRNAFasta);
    if (retVal != 0)
    {
        return 1;
    }

    // Read mRNA fasta file
    retVal = mMRNASeqs.read_fasta(mMRNAFasta);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

//
// MR3Core methods
//
template <class TRNAString, int SEEDLEN>
void MR3Core<TRNAString, SEEDLEN>::init_from_args(MR3Options& opts)
{
    mOutputAlign = opts.mOutputAlign;
    mOFileSite = opts.mOFileSite;
    mOFileTotal = opts.mOFileTotal;
    mMinSeedLen = opts.mMinSeedLen;
    mMaxSeedLen = opts.mMaxSeedLen;
    mMinAlignScore = opts.mMinAlignScore;
    mMaxEnergy = opts.mMaxEnergy;

    resize(mSeedDef, 6);
    mSeedDef[0] = 'Y';
    mSeedDef[1] = 'Y';
    mSeedDef[2] = 'Y';
    if (opts.mMinSeedLen == 7)
    {
        mSeedDef[0] = 'N';
    }
    else if (opts.mMinSeedLen == 8)
    {
        mSeedDef[0] = 'N';
        mSeedDef[1] = 'N';
    }

    if (opts.mMaxSeedLen == 7)
    {
        mSeedDef[2] = 'N';
    }
    else if (opts.mMaxSeedLen == 6)
    {
        mSeedDef[2] = 'N';
        mSeedDef[1] = 'N';
    }
    mSeedDef[3] = opts.mAllowGUWobble;
    mSeedDef[4] = opts.mAllowMismatch;
    mSeedDef[5] = opts.mAllowBT;

    mSiteScores.set_min_align_score(mMinAlignScore);
    mSiteScores.set_max_energy(mMaxEnergy);

}

template <class TRNAString, int SEEDLEN>
int MR3Core<TRNAString, SEEDLEN>::open_output_file()
{
    // Open output file 1
    mOFile1.open(toCString(mOFileSite), std::ofstream::out);
    if (!mOFile1.good())
    {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileSite) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Open output file 2
    mOFile2.open(toCString(mOFileTotal), std::ofstream::out);
    if (!mOFile2.good())
    {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileTotal) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return 0;
}

template <class TRNAString, int SEEDLEN>
int MR3Core<TRNAString, SEEDLEN>::calculate_all_scores()
{
    int retVal;

    for (unsigned i = 0; i < length(mMiRNASeqs); ++i)
    {

#if SEQAN_ENABLE_DEBUG
        clock_t startTime = clock();
#endif

        retVal = calculate_mirna_scores(i);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Score calculation failed for ";
            std::cerr << toCString((seqan::CharString)mMiRNAIds[i]) << "." << std::endl;
            return 1;
        }

#if SEQAN_ENABLE_DEBUG
        std::cout << toCString((seqan::CharString)mMiRNAIds[i]) << ": ";
        std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;
#endif

    }

    return 0;
}

template <class TRNAString, int SEEDLEN>
int MR3Core<TRNAString, SEEDLEN>::calculate_mirna_scores(unsigned pIdx)
{
    int retVal;
    TRNAString miRNASeq = mMiRNASeqs[pIdx];

    // Search seed sites
    if (mExecSearchSeedSites)
    {
        retVal = mSeedSites.find_seed_sites(miRNASeq, mSeedDef);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Calculate alignment and energy scores
    if (mExecCalSiteScore)
    {
        retVal = mSiteScores.calc_scores(mSeedSites, miRNASeq, mMRNASeqs);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Calculate site scores failed." << std::endl;
            return 1;
        }
    }

    // Filter overlapped sites
    if (OVERLAP_LEN != 0 && mExecFilterOverlap)
    {
        retVal = mOverlappedSites.filter_overlapped_sites_by_scores(mSeedSites, mSiteScores, OVERLAP_LEN);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Check overlapped sites failed." << std::endl;
            return 1;
        }
    }


    // Sort target sites
    if (mExecSortSites)
    {
        retVal = mSortedSites.generate_sorted_mrna_pos(mSeedSites);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Sort target sites failed." << std::endl;
            return 1;
        }
    }

    // Summarize alignment and energy scores
    if (mExecSumScores)
    {
        const seqan::String<unsigned>& sortedPos = mSortedSites.get_sorted_mrna_pos();
        retVal = mTotalScores.calc_scores(mSeedSites, mSiteScores, sortedPos);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Calculate total scores failed." << std::endl;
            return 1;
        }
    }

    // Write site scores
    if (mOutputSiteScore)
    {
        retVal = write_site_score(mMiRNAIds[pIdx]);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Could not write site scores." << std::endl;
            return 1;
        }
    }

    // Write total scores
    if (mOutputTotalScore)
    {
        retVal = write_total_score(mMiRNAIds[pIdx]);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Could not write total ddG scores." << std::endl;
            return 1;
        }
    }

    // Write alignments
    if (mOutputAlign)
    {
        retVal = write_alignment(mMiRNAIds[pIdx]);
        if (retVal != 0)
        {
            std::cerr << "ERROR: Could not write alignments." << std::endl;
            return 1;
        }
    }

    mSeedSites.clear_pos();
    mOverlappedSites.clear_cluster();
    mSiteScores.clear_scores();
    mSortedSites.clear_site_pos();
    mTotalScores.clear_scores();

    return 0;
}

template <class TRNAString, int SEEDLEN>
int MR3Core<TRNAString, SEEDLEN>::write_site_score(seqan::CharString const &pMiRNAId)
{
    const seqan::String<unsigned>& mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned>& sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString>& seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned>& sortedPos = mSortedSites.get_sorted_mrna_pos();
    seqan::CharString seedType;
    int seedStart = 0;
    int posIdx;
    float scoreAlign;
    float scoreEn;

    for (unsigned i = 0; i < length(sortedPos); ++i)
    {
        posIdx = sortedPos[i];
        if (!mSeedSites.mEffectiveSites[posIdx])
        {
            continue;
        }

        seedStart = sitePos[posIdx];
        scoreAlign =  mSiteScores.get_align_score(posIdx);
        scoreAlign = roundf(scoreAlign * 100.0f) / 100.0f;
        scoreEn =  mSiteScores.get_energy_score(posIdx);
        scoreEn = roundf(scoreEn * 100.0f) / 100.0f;

        mOFile1 << toCString(pMiRNAId) << "\t";
        mOFile1 << toCString((seqan::CharString)mMRNAIds[mRNAPos[posIdx]]) << "\t";
        mOFile1 << seedStart + 1  << "\t";
        mOFile1 << seedStart + 1 + INDEXED_SEQ_LEN << "\t";
        mOFile1 << toCString((seqan::CharString)seedTypes[posIdx]) << "\t";
        mOFile1 << scoreAlign << "\t";
        mOFile1 << scoreEn << "\t";
        mOFile1 << std::endl;
    }

    return 0;
}

template <class TRNAString, int SEEDLEN>
int MR3Core<TRNAString, SEEDLEN>::write_total_score(seqan::CharString const &pMiRNAId)
{
    const seqan::String<float>& totalAlignScores = mTotalScores.get_align_scores();
    const seqan::String<float>& totalEnScores = mTotalScores.get_energy_scores();
    const seqan::String<int>& mRNAPos = mTotalScores.get_mrna_pos();
    const seqan::String<int>& siteNum = mTotalScores.get_site_num();

    for (unsigned i = 0; i < length(mRNAPos); ++i)
    {

        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString)mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile2 << totalAlignScores[i] << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << totalEnScores[i] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

template <class TRNAString, int SEEDLEN>
int MR3Core<TRNAString, SEEDLEN>::write_alignment(seqan::CharString const &pMiRNAId)
{
    const seqan::String<unsigned>& mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned>& sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString>& seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned>& sortedPos = mSortedSites.get_sorted_mrna_pos();
    seqan::CharString seedType;
    int seedStart = 0;
    int posIdx;
    int count = 0;
    float align_score;
    float energy_score;

    for (unsigned i = 0; i < length(sortedPos); ++i)
    {
        posIdx = sortedPos[i];

        if (!mSeedSites.mEffectiveSites[posIdx])
        {
            continue;
        }

        seedStart = sitePos[posIdx];
        align_score = mSiteScores.get_align_score(posIdx);
        align_score = roundf(align_score * 100.0f) / 100.0f;
        energy_score = mSiteScores.get_energy_score(posIdx);
        energy_score = roundf(energy_score * 100.0f) / 100.0f;

        std::cout << "### " << count+1 << ": " << toCString(pMiRNAId) <<" ###" << std::endl;
        mSiteScores.print_alignment(posIdx);
        std::cout << "  miRNA:           " << toCString(pMiRNAId) << std::endl;
        std::cout << "  mRNA:            " << toCString((seqan::CharString)mMRNAIds[mRNAPos[posIdx]]) << std::endl;
        std::cout << "  seed type:       " << toCString((seqan::CharString)seedTypes[posIdx]) << std::endl;
        std::cout << "  position(start): " << seedStart + 1 << std::endl;
        std::cout << "  position(end):   " << seedStart + 1 + INDEXED_SEQ_LEN << std::endl;
        std::cout << "  alignment score: " << align_score << std::endl;
        std::cout << "  energy score:    " << energy_score << std::endl;

        std::cout << std::endl;

        ++count;

    }

    return 0;
}

// Explicit template instantiation
template class MR3CoreInput<TRNATYPE>;
template class MR3Core<TRNATYPE>;

} // namespace mr3as
