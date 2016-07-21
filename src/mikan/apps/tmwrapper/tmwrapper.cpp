#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <map>
#include <unistd.h>
#include <sys/param.h>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

struct TMWOptions
{
    unsigned max_pair;
    unsigned max_seq_len;
    CharString tm_path;
    CharString tm_pos_path;
    CharString tm_exe;
    CharString data_path;
    CharString ofile_mirna_fa;
    CharString ofile_mrna_fa;
    CharString ofile_mirna_mrna;
    CharString mirna_fa;
    CharString mrna_fa;
    CharString ofile;
    bool noexec;
    bool exec_pos_version;

    TMWOptions() : max_pair(900), max_seq_len(4999),
            tm_path("/scratch/tsa103/workspace/uib.mirna_target_eval/mirna_target_eval/programs/TargetMiner/score/"),
            tm_pos_path("/scratch/tsa103/workspace/uib.mirna_target_eval/mirna_target_eval/programs/TargetMiner/pos/"),
            tm_exe("TargetMiner"),
            data_path("/scratch/tsa103/data/seqan/mirna_target/targetminer/"),
            ofile_mirna_fa("tmp_mirna.fa"),
            ofile_mrna_fa("tmp_mrna.fa"),
            ofile_mirna_mrna("tmp_mirna_mrna.txt"),
            noexec (false), exec_pos_version (false)
    {}
};

ArgumentParser::ParseResult parseCommandLine(
        TMWOptions &options,
        int const argc,
        char const **argv)
{
    // Setup ArgumentParser
    ArgumentParser parser("tmwrapper");

    // Define Options
    addSection(parser, "TargetMiner Wrapper Options");
    addOption(parser, ArgParseOption("n", "noexec", "Do not execute TargetMinder but create necessary files."));
    addOption(parser, ArgParseOption("p", "outpos", "Output positions instead of scores."));

    // Set short description, version, and date
    setShortDescription(parser, "Wrapper program for TargetMiner");
    setVersion(parser, "1.0");
    setDate(parser, "January 2014");

    // Define usage line and long description
    addUsageLine(parser,
            "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" \"\\fIOUTPUT FILE\\fP\"");
    addDescription(parser,
            "This program creates input files necessary for TargetMiner and executes TargetMiner with them.");

    // Define Arguments
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
            "\\fBtmwrapper\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP \\fIoutput_file\\fP",
            "create input files from \\fImiRNAs\\fP and \\fImRNA\\fP, call TargetMinder with them, "
            "and append the results to \\fIoutput\\fP.");

    // Parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        return res;
    }

    // Extract arguments
    getArgumentValue(options.mirna_fa, parser, 0);
    getArgumentValue(options.mrna_fa, parser, 1);
    getArgumentValue(options.ofile, parser, 2);
    std::fstream mirna_fa(toCString(options.mirna_fa), std::ios::binary | std::ios::in);
    if (!mirna_fa.good())
    {
        std::cerr << "ERROR: Could not open input file " << options.mirna_fa << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    std::fstream mrna_fa(toCString(options.mrna_fa), std::ios::binary | std::ios::in);
    if (!mrna_fa.good())
    {
        std::cerr << "ERROR: Could not open input file " << options.mrna_fa << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    std::fstream ofile(toCString(options.ofile), std::ios::binary | std::ios::out);
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << options.ofile << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Extract options
    options.noexec = isSet(parser, "noexec");
    options.exec_pos_version = isSet(parser, "outpos");

    return ArgumentParser::PARSE_OK;
}

int read_fasta(
        CharString const &fasta,
        StringSet<CharString> &ids,
        StringSet<CharString> &seqs,
        bool pIsMRNA)
{
    std::multimap<unsigned, unsigned> mrna_map;
    std::multimap<unsigned, unsigned>::iterator it;
    unsigned i = 0;
    unsigned j = 0;
    CharString id;
    CharString seq;
    StringSet<CharString> tmp_ids;
    StringSet<CharString> tmp_seqs;

    SequenceStream seqStream(toCString(fasta));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file!" << std::endl;
        return 1;
    }

    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from " << fasta << "!" << std::endl;
            return 1;
        }

        if (pIsMRNA)
        {
            toLower(seq);
        }
        else
        {
            toUpper(seq);
        }

        for (unsigned i = 0; i < length(seq); ++i)
        {
            if (pIsMRNA)
            {
                if (seq[i] == 'u')
                {
                    seq[i] = 't';
                }
            }
            else
            {
                if (seq[i] == 'T')
                {
                    seq[i] = 'U';
                }
            }
        }

        if (!pIsMRNA)
        {
            appendValue(seqs, seq);
            appendValue(ids, id);
        }
        else
        {
            appendValue(tmp_seqs, seq);
            appendValue(tmp_ids, id);
        }
    }

    if (!pIsMRNA)
    {
        return 0;
    }

    mrna_map.clear();
    for (unsigned i = 0; i < length(tmp_seqs); ++i)
    {
        mrna_map.insert(std::pair<unsigned, unsigned>(length(tmp_seqs[i]), i));
    }

    resize(seqs, length(tmp_seqs));
    resize(ids, length(tmp_ids));
    i = 0;
    j = 0;
    for (it = mrna_map.begin(); it != mrna_map.end(); ++it)
    {
        j = (*it).second;
//        std::cout << i << "," << j << "," << length(tmp_seqs[j]) << std::endl;
        seqs[i] = tmp_seqs[j];
        ids[i] = tmp_ids[j];
        ++i;
    }

    return 0;
}

int create_mirna_fa(
        TMWOptions& pOptions,
        StringSet<CharString> const &pMiRNAIds,
        StringSet<CharString> const &pMiRNAs,
        unsigned pI)
{

    CharString ofname = pOptions.data_path;
    append(ofname, pOptions.ofile_mirna_fa);

    std::fstream ofile(toCString(ofname), std::ios::binary | std::ios::out);
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << ofname << std::endl;
        return 1;
    }

    ofile << ">" << pMiRNAIds[pI] << "\t" << pMiRNAIds[pI] << std::endl;
    ofile << pMiRNAs[pI] << std::endl;
    ofile << ">" << std::endl;

    ofile.close();

    return 0;
}

int create_mrna_fa(
        TMWOptions& pOptions,
        StringSet<CharString> const &pMNRAIds,
        StringSet<CharString> &pMRNAs,
        unsigned pStartJ,
        unsigned pEndJ)
{
    CharString ofname = pOptions.data_path;
    append(ofname, pOptions.ofile_mrna_fa);

    std::fstream ofile(toCString(ofname), std::ios::binary | std::ios::out);
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << ofname << std::endl;
        return 1;
    }

    for (unsigned i = pStartJ; i < pEndJ + 1; ++i)
    {
        ofile << ">hg18 " << pMNRAIds[i] << std::endl;
        ofile << pMRNAs[i] << std::endl;
    }
    ofile << ">" << std::endl;

    ofile.close();

    return 0;
}

int create_mirna_mrna_list(
        TMWOptions& pOptions,
        StringSet<CharString> const &pMiRNAIds,
        StringSet<CharString> const &pMNRAIds,
        StringSet<CharString>& pList,
        unsigned pI,
        unsigned pStartJ,
        unsigned pEndJ)
{
    unsigned idx;
    CharString ofname = pOptions.data_path;
    append(ofname, pOptions.ofile_mirna_mrna);

    std::fstream ofile(toCString(ofname), std::ios::binary | std::ios::out);
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << ofname << std::endl;
        return 1;
    }

    clear(pList);
    resize(pList, pEndJ + 1 - pStartJ);
    idx = 0;
    for (unsigned j = pStartJ; j < pEndJ + 1; ++j)
    {
        ofile << ">" << pMiRNAIds[pI] << "\t" << pMNRAIds[j] << std::endl;
        append(pList[idx], pMiRNAIds[pI]);
        append(pList[idx], "\t");
        append(pList[idx], pMNRAIds[j]);
        append(pList[idx], "\t");
        ++idx;
    }
    ofile << ">" << std::endl;

    ofile.close();

    return 0;
}

int write_pos_output(TMWOptions& pOptions, std::fstream &pOfile)
{
    std::string line;
    CharString ifname = pOptions.tm_pos_path;
    append(ifname, "TargetMiner_Prediction.html");
    std::string start_str ("Predicted miRNA Targets");
    std::size_t found;

    std::ifstream ifile(toCString(ifname), std::ios::in);
    if (!ifile.good())
    {
        std::cerr << "ERROR: Could not open result file " << ifname << std::endl;
        return 1;
    }

    while (getline(ifile, line))
    {
        found = line.find(start_str);
        if (found != std::string::npos)
        {
            break;
        }
    }

    pOfile << "#pos_start" << std::endl;
    while (getline(ifile, line))
    {
        pOfile << line << std::endl;
    }
    pOfile << "#pos_end" << std::endl;

    ifile.close();

    return 0;
}

int rw_output(TMWOptions& pOptions, StringSet<CharString>& pList, std::fstream &pOfile)
{
    int label;
    double dec_val;
    unsigned idx;

    CharString ifname = pOptions.tm_path;
    append(ifname, "output");

    std::ifstream ifile(toCString(ifname), std::ios::in);
    if (!ifile.good())
    {
        std::cerr << "ERROR: Could not open result file " << ifname << std::endl;
        return 1;
    }

    idx = 0;
    while (ifile >> label >> dec_val)
    {
        pOfile << pList[idx] << label << "\t" << dec_val << std::endl;
        ++idx;
    }

    ifile.close();

    return 0;
}

int write_output(TMWOptions&, StringSet<CharString>& pList, std::fstream &pOfile)
{
    for (unsigned i = 0; i < length(pList); ++i)
    {
        pOfile << pList[i] << "1" << "\t" << "NA" << std::endl;
    }

    return 0;
}

int call_targetminer(TMWOptions& pOptions, StringSet<CharString>& pList, std::fstream &pOfile)
{
    int ret;

    if (pOptions.noexec)
    {
        ret = write_output(pOptions, pList, pOfile);
        return 0;
    }

    CharString ifname_mirna = pOptions.data_path;
    CharString ifname_mrna = pOptions.data_path;
    CharString ifname_list = pOptions.data_path;

    append(ifname_mirna, pOptions.ofile_mirna_fa);
    append(ifname_mrna, pOptions.ofile_mrna_fa);
    append(ifname_list, pOptions.ofile_mirna_mrna);

    CharString exe_str;
    if (pOptions.exec_pos_version)
    {
        exe_str = pOptions.tm_pos_path;
    }
    else
    {
        exe_str = pOptions.tm_path;
    }
    append(exe_str, pOptions.tm_exe);
    append(exe_str, " -test ");
    append(exe_str, ifname_list);
    append(exe_str, " -3utr ");
    append(exe_str, ifname_mrna);
    append(exe_str, " -mir ");
    append(exe_str, ifname_mirna);

//    std::cout << exe_str << std::endl;
    ret = std::system(toCString(exe_str));
    if (!pOptions.exec_pos_version)
    {
        if (ret)
        {
    //        std::cerr << "ERROR: Failed to call TargetMiner " << std::endl;
            ret = write_output(pOptions, pList, pOfile);
        }
        else
        {
            ret = rw_output(pOptions, pList, pOfile);
        }
    }
    else
    {
        ret = write_pos_output(pOptions, pOfile);
    }

    if (ret)
    {
        std::cerr << "ERROR: Failed to write output." << std::endl;
        return ret;
    }

    return 0;
}

int create_files_and_call_prog(
        TMWOptions& pOptions,
        StringSet<CharString> const &pMiRNAIds,
        StringSet<CharString> const &pMiRNAs,
        StringSet<CharString> const &pMNRAIds,
        StringSet<CharString> &pMRNAs,
        unsigned pI,
        unsigned pStartJ,
        unsigned pEndJ,
        std::fstream &pOfile)
{
    int ret;
    StringSet<CharString> mirna_mrnas;

    ret = create_mirna_fa(pOptions, pMiRNAIds, pMiRNAs, pI);
    if (ret)
    {
        return ret;
    }

    ret = create_mrna_fa(pOptions, pMNRAIds, pMRNAs, pStartJ, pEndJ);
    if (ret)
    {
        return ret;
    }

    ret = create_mirna_mrna_list(pOptions, pMiRNAIds, pMNRAIds, mirna_mrnas, pI, pStartJ, pEndJ);
    if (ret)
    {
        return ret;
    }

    ret = call_targetminer(pOptions, mirna_mrnas, pOfile);
    if (ret)
    {
        return ret;
    }

    return 0;
}

int exec_tm(
        TMWOptions& pOptions,
        StringSet<CharString> const &pMiRNAIds,
        StringSet<CharString> const &pMiRNAs,
        StringSet<CharString> const &pMNRAIds,
        StringSet<CharString> &pMRNAs)
{
    unsigned startJ, endJ;
    unsigned pairCount;
    bool execNext;
    int ret;

    char *pathname = get_current_dir_name();
    CharString path_str;
    if (pOptions.exec_pos_version)
    {
        path_str = pOptions.tm_pos_path;
    }
    else
    {
        path_str = pOptions.tm_path;
    }
    if (chdir(toCString(path_str)) == -1)
    {
        std::cerr << "ERROR: Failed to change directory " << std::endl;
        return 1;
    }

    std::fstream ofile(toCString(pOptions.ofile), std::ios::binary | std::ios::out);
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << pOptions.ofile << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(pMiRNAIds); ++i)
    {
        startJ = 0;
        endJ = 0;
        pairCount = 0;
        execNext = false;
        for (unsigned j = 0; j < length(pMNRAIds); ++j)
        {
            if (execNext)
            {
//                std::cout << i << ", " << startJ << " - " << endJ << std::endl;
                ret = create_files_and_call_prog(pOptions, pMiRNAIds, pMiRNAs, pMNRAIds, pMRNAs, i, startJ, endJ,
                        ofile);
                if (ret)
                {
                    return ret;
                }
                startJ = j;

                execNext = false;
                pairCount = 0;
            }

            ++pairCount;
            if(pairCount == pOptions.max_pair || length(pMRNAs[j]) > pOptions.max_seq_len)
            {
                endJ = j;
                execNext = true;
            }
        }

        endJ = length(pMNRAIds) - 1;
//        std::cout << i << ", " << startJ << " - " << endJ << std::endl;
        ret = create_files_and_call_prog(pOptions, pMiRNAIds, pMiRNAs, pMNRAIds, pMRNAs, i, startJ, endJ, ofile);
        if (ret)
        {
            return ret;
        }
    }

    if (chdir(pathname) == -1)
    {
        std::cerr << "ERROR: Failed to change directory " << std::endl;
        return 1;
    }

    ofile.close();

    return 0;
}

int main(int argc, char const ** argv)
{
    // Declare variables
    StringSet<CharString> mirna_ids;
    StringSet<CharString> mrna_ids;
    StringSet<CharString> mirna_seqs;
    StringSet<CharString> mrna_seqs;
    int ret;

    // Parse the command line.
    TMWOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read fasta files
    int fread_res;
    fread_res = read_fasta(options.mirna_fa, mirna_ids, mirna_seqs, false);
    if (fread_res != 0)
    {
        return fread_res;
    }
    fread_res = read_fasta(options.mrna_fa, mrna_ids, mrna_seqs, true);
    if (fread_res != 0)
    {
        return fread_res;
    }

    ret = exec_tm(options, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs);
    if (ret)
    {
        return ret;
    }

    return 0;
}
