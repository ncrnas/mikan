#ifndef MKE_CONFIG_HPP_
#define MKE_CONFIG_HPP_

#include <seqan/sequence.h>
#include "inih/INIReader.h"
#include "mke_enum.hpp"           // ToolIdx, SiteIdx, RNAIdx

namespace mkens {

//
// Tool config interface
//
class MKEConfig {
public:
    static const int CommentPos = 26;

    // Define method
    MKEConfig() :
            mTrue("true"),
            mFalse("false") {

        //
        // Basic parameters - tool level
        //
        // Tool keys
        mToolKeys.resize(ToolIdx::Count);
        mToolKeys[ToolIdx::MR] = "mr";
        mToolKeys[ToolIdx::PT] = "pt";
        mToolKeys[ToolIdx::RH] = "rh";
        mToolKeys[ToolIdx::TM] = "tm";
        mToolKeys[ToolIdx::TS] = "ts";
        mToolKeys[ToolIdx::SV] = "sv";

        // Tool descriptions
        mToolDesc.resize(ToolIdx::Count);
        mToolDesc[ToolIdx::MR] = "miRanda 3";
        mToolDesc[ToolIdx::PT] = "PITA";
        mToolDesc[ToolIdx::RH] = "RNAhybrid";
        mToolDesc[ToolIdx::TM] = "TargetMiner";
        mToolDesc[ToolIdx::TS] = "TargetScan 5";
        mToolDesc[ToolIdx::SV] = "Two-step SVM";

        // Tool flags
        mToolDefFlg.resize(ToolIdx::Count);
        mToolDefFlg[ToolIdx::MR] = true;
        mToolDefFlg[ToolIdx::PT] = true;
        mToolDefFlg[ToolIdx::RH] = true;
        mToolDefFlg[ToolIdx::TM] = true;
        mToolDefFlg[ToolIdx::TS] = true;
        mToolDefFlg[ToolIdx::SV] = true;

        //
        // Basic parameters - site level
        //
        // Score keys
        mSiteKeys.resize(SiteIdx::Count);
        mSiteKeys[SiteIdx::MRAlg] = "mr:site:alg";
        mSiteKeys[SiteIdx::MREng] = "mr:site:eng";
        mSiteKeys[SiteIdx::PTDdg] = "pt:site:ddg";
        mSiteKeys[SiteIdx::PTDpx] = "pt:site:dpx";
        mSiteKeys[SiteIdx::PTOpn] = "pt:site:opn";
        mSiteKeys[SiteIdx::RHMfe] = "rh:site:mfe";
        mSiteKeys[SiteIdx::RHNrm] = "rh:site:nrm";
        mSiteKeys[SiteIdx::TSCtx] = "ts:site:ctx";
        mSiteKeys[SiteIdx::SVSvm] = "sv:site:svm";

        // Score descriptions
        mSiteDesc.resize(SiteIdx::Count);
        mSiteDesc[SiteIdx::MRAlg] = "miRanda alignment score";
        mSiteDesc[SiteIdx::MREng] = "miRanda energy score";
        mSiteDesc[SiteIdx::PTDdg] = "PITA ddG";
        mSiteDesc[SiteIdx::PTDpx] = "PITA dDuplex";
        mSiteDesc[SiteIdx::PTOpn] = "PITA dOpen";
        mSiteDesc[SiteIdx::RHMfe] = "RNAhybrid MFE";
        mSiteDesc[SiteIdx::RHNrm] = "RNAhybrid MFE normalized by length";
        mSiteDesc[SiteIdx::TSCtx] = "TargetScan context score";
        mSiteDesc[SiteIdx::SVSvm] = "Two-step SVM score";

        // Score flags
        mSiteDefFlg.resize(SiteIdx::Count);
        mSiteDefFlg[SiteIdx::MRAlg] = true;
        mSiteDefFlg[SiteIdx::MREng] = true;
        mSiteDefFlg[SiteIdx::PTDdg] = true;
        mSiteDefFlg[SiteIdx::PTDpx] = true;
        mSiteDefFlg[SiteIdx::PTOpn] = true;
        mSiteDefFlg[SiteIdx::RHMfe] = true;
        mSiteDefFlg[SiteIdx::RHNrm] = true;
        mSiteDefFlg[SiteIdx::TSCtx] = true;
        mSiteDefFlg[SiteIdx::SVSvm] = true;

        //
        // Basic parameters - RNA level
        //
        // Score keys
        mRNAKeys.resize(RNAIdx::Count);
        mRNAKeys[RNAIdx::MRAlg] = "mr:rna:alg";
        mRNAKeys[RNAIdx::MREng] = "mr:rna:eng";
        mRNAKeys[RNAIdx::PTDdg] = "pt:rna:ddg";
        mRNAKeys[RNAIdx::RHMfe] = "rh:rna:mfe";
        mRNAKeys[RNAIdx::RHNrm] = "rh:rna:nrm";
        mRNAKeys[RNAIdx::TMSvm] = "tm:rna:svm";
        mRNAKeys[RNAIdx::TSCtx] = "ts:rna:ctx";
        mRNAKeys[RNAIdx::SVSvm] = "sv:rna:svm";

        // Score descriptions
        mRNADesc.resize(RNAIdx::Count);
        mRNADesc[RNAIdx::MRAlg] = "max + log(total) miRanda alignment score";
        mRNADesc[RNAIdx::MREng] = "min + log(total) miRanda energy score";
        mRNADesc[RNAIdx::PTDdg] = "PITA ddG";
        mRNADesc[RNAIdx::RHMfe] = "min + log(total) RNAhybrid MFE";
        mRNADesc[RNAIdx::RHNrm] = "max + log(total) RNAhybrid MFE normalized by length";
        mRNADesc[RNAIdx::TMSvm] = "TargetMiner SVM";
        mRNADesc[RNAIdx::TSCtx] = "TargetScan context score";
        mRNADesc[RNAIdx::SVSvm] = "Two-step SVM score";

        // Score flags
        mRNADefFlg.resize(RNAIdx::Count);
        mRNADefFlg[RNAIdx::MRAlg] = true;
        mRNADefFlg[RNAIdx::MREng] = true;
        mRNADefFlg[RNAIdx::PTDdg] = true;
        mRNADefFlg[RNAIdx::RHMfe] = true;
        mRNADefFlg[RNAIdx::RHNrm] = true;
        mRNADefFlg[RNAIdx::TMSvm] = true;
        mRNADefFlg[RNAIdx::TSCtx] = true;
        mRNADefFlg[RNAIdx::SVSvm] = true;

        //
        // Default weight values
        //
        // Site level
        mSiteDefWeight.resize(SiteIdx::Count);
        mSiteDefWeight[SiteIdx::MRAlg] = 6;
        mSiteDefWeight[SiteIdx::MREng] = 6;
        mSiteDefWeight[SiteIdx::PTDdg] = 8;
        mSiteDefWeight[SiteIdx::PTDpx] = 2;
        mSiteDefWeight[SiteIdx::PTOpn] = 2;
        mSiteDefWeight[SiteIdx::RHMfe] = 9;
        mSiteDefWeight[SiteIdx::RHNrm] = 3;
        mSiteDefWeight[SiteIdx::TSCtx] = 12;
        mSiteDefWeight[SiteIdx::SVSvm] = 12;

        // RNA level
        mRNADefWeight.resize(RNAIdx::Count);
        mRNADefWeight[RNAIdx::MRAlg] = 6;
        mRNADefWeight[RNAIdx::MREng] = 6;
        mRNADefWeight[RNAIdx::PTDdg] = 12;
        mRNADefWeight[RNAIdx::RHMfe] = 9;
        mRNADefWeight[RNAIdx::RHNrm] = 3;
        mRNADefWeight[RNAIdx::TMSvm] = 12;
        mRNADefWeight[RNAIdx::TSCtx] = 12;
        mRNADefWeight[RNAIdx::SVSvm] = 12;

        //
        // Default normalization parameters - site level
        //
        mSiteDefLower.resize(SiteIdx::Count);
        mSiteDefUpper.resize(SiteIdx::Count);
        mSiteDefReverse.resize(SiteIdx::Count);

        // miRanda alignment score
        mSiteDefLower[SiteIdx::MRAlg] = 140;
        mSiteDefUpper[SiteIdx::MRAlg] = 200;
        mSiteDefReverse[SiteIdx::MRAlg] = false;

        // miRanda energy score
        mSiteDefLower[SiteIdx::MREng] = -55;
        mSiteDefUpper[SiteIdx::MREng] = 0;
        mSiteDefReverse[SiteIdx::MREng] = true;

        // PITA ddG
        mSiteDefLower[SiteIdx::PTDdg] = -42;
        mSiteDefUpper[SiteIdx::PTDdg] = 36;
        mSiteDefReverse[SiteIdx::PTDdg] = true;

        // PITA dDuplex
        mSiteDefLower[SiteIdx::PTDpx] = -48;
        mSiteDefUpper[SiteIdx::PTDpx] = 0;
        mSiteDefReverse[SiteIdx::PTDpx] = true;

        // PITA dOpen
        mSiteDefLower[SiteIdx::PTOpn] = -47;
        mSiteDefUpper[SiteIdx::PTOpn] = 0;
        mSiteDefReverse[SiteIdx::PTOpn] = true;

        // RNAhybrid MFE
        mSiteDefLower[SiteIdx::RHMfe] = -42;
        mSiteDefUpper[SiteIdx::RHMfe] = 0;
        mSiteDefReverse[SiteIdx::RHMfe] = true;

        // RNAhybrid MFE normalized by length
        mSiteDefLower[SiteIdx::RHNrm] = 0;
        mSiteDefUpper[SiteIdx::RHNrm] = 8;
        mSiteDefReverse[SiteIdx::RHNrm] = false;

        // TargetScan context score
        mSiteDefLower[SiteIdx::TSCtx] = -0.64f;
        mSiteDefUpper[SiteIdx::TSCtx] = 0.2;
        mSiteDefReverse[SiteIdx::TSCtx] = false;

        // Two-step SVM score"
        mSiteDefLower[SiteIdx::SVSvm] = -2.54f;
        mSiteDefUpper[SiteIdx::SVSvm] = 3.06f;
        mSiteDefReverse[SiteIdx::SVSvm] = false;

        //
        // Default normalization parameters - RNA level
        //
        // RNA score default parameters for normalization
        mRNADefLower.resize(RNAIdx::Count);
        mRNADefUpper.resize(RNAIdx::Count);
        mRNADefReverse.resize(RNAIdx::Count);

        // max + log(total) miRanda alignment score
        mRNADefLower[RNAIdx::MRAlg] = 140;
        mRNADefUpper[RNAIdx::MRAlg] = 210;
        mRNADefReverse[RNAIdx::MRAlg] = false;

        // tmin + log(total) miRanda energy score
        mRNADefLower[RNAIdx::MREng] = -62;
        mRNADefUpper[RNAIdx::MREng] = 1;
        mRNADefReverse[RNAIdx::MREng] = true;

        // PITA ddG
        mRNADefLower[RNAIdx::PTDdg] = -48;
        mRNADefUpper[RNAIdx::PTDdg] = 30;
        mRNADefReverse[RNAIdx::PTDdg] = true;

        // min + log(total) RNAhybrid MFE
        mRNADefLower[RNAIdx::RHMfe] = -50;
        mRNADefUpper[RNAIdx::RHMfe] = 0;
        mRNADefReverse[RNAIdx::RHMfe] = true;

        // max + log(total) RNAhybrid MFE normalized by length
        mRNADefLower[RNAIdx::RHNrm] = 0;
        mRNADefUpper[RNAIdx::RHNrm] = 13;
        mRNADefReverse[RNAIdx::RHNrm] = false;

        // TargetMiner SVM
        mRNADefLower[RNAIdx::TMSvm] = -1.4f;
        mRNADefUpper[RNAIdx::TMSvm] = 1;
        mRNADefReverse[RNAIdx::TMSvm] = false;

        // TargetScan context score
        mRNADefLower[RNAIdx::TSCtx] = -13;
        mRNADefUpper[RNAIdx::TSCtx] = 1.8f;
        mRNADefReverse[RNAIdx::TSCtx] = false;

        // Two-step SVM score
        mRNADefLower[RNAIdx::SVSvm] = -8;
        mRNADefUpper[RNAIdx::SVSvm] = 20;
        mRNADefReverse[RNAIdx::SVSvm] = false;

        init_config();
    }

    std::string& get_tool_key(unsigned pIdx) {
        return mToolKeys[pIdx];
    }

    bool get_tool_flag(unsigned pIdx) const {
        return get_bool_from_map(mToolFlg, mToolKeys[pIdx]);
    }

    bool get_site_flag(std::string &pKey) const {
        return get_bool_from_map(mSiteFlg, pKey);
    }

    bool get_rna_flag(std::string &pKey) const {
        return get_bool_from_map(mRNAFlg, pKey);
    }

    float get_site_lower(std::string &pKey) const {
        return get_float_from_map(mSiteLower, pKey);
    }

    float get_site_upper(std::string &pKey) const {
        return get_float_from_map(mSiteUpper, pKey);
    }

    bool get_site_reverse(std::string &pKey) const {
        return get_bool_from_map(mSiteReverse, pKey);
    }

    float get_site_weight(std::string &pKey) const {
        return get_float_from_map(mSiteWeight, pKey);
    }

    float get_rna_lower(std::string &pKey) const {
        return get_float_from_map(mRNALower, pKey);
    }

    float get_rna_upper(std::string &pKey) const {
        return get_float_from_map(mRNAUpper, pKey);
    }

    float get_rna_reverse(std::string &pKey) const {
        return get_bool_from_map(mRNAReverse, pKey);
    }

    float get_rna_weight(std::string &pKey) const {
        return get_float_from_map(mRNAWeight, pKey);
    }

    // Method prototypes
    void init_config();

    void print_config();

    int parse_config(std::string &pConfFile);

private:
    typedef std::vector<std::string>::iterator TKeyIt;

    void init_flags();

    void init_weights();

    void init_tool_config();

    std::string &get_bool_str(bool pFlag);

    void print_flags();

    void print_weights();

    void print_tool_config();

    int parse_flags(INIReader &pReader);

    int parse_weights(INIReader &pReader);

    int parse_tool_config(INIReader &pReader);

    // Flags
    std::map<std::string, bool> mToolFlg;
    std::map<std::string, bool> mSiteFlg;
    std::map<std::string, bool> mRNAFlg;

    // Weight values
    std::map<std::string, float> mSiteWeight;
    std::map<std::string, float> mRNAWeight;

    // Normalization parameters - site level
    std::map<std::string, float> mSiteLower;
    std::map<std::string, float> mSiteUpper;
    std::map<std::string, bool> mSiteReverse;

    // Normalization parameters - RNA level
    std::map<std::string, float> mRNALower;
    std::map<std::string, float> mRNAUpper;
    std::map<std::string, bool> mRNAReverse;

    // Basic parameters - tool level
    std::vector<std::string> mToolKeys;
    std::vector<std::string> mToolDesc;
    std::vector<bool> mToolDefFlg;

    // Basic parameters - site level
    std::vector<std::string> mSiteKeys;
    std::vector<std::string> mSiteDesc;
    std::vector<bool> mSiteDefFlg;

    // Basic parameters - RNA level
    std::vector<std::string> mRNAKeys;
    std::vector<std::string> mRNADesc;
    std::vector<bool> mRNADefFlg;

    // Default weight values
    std::vector<float> mSiteDefWeight;
    std::vector<float> mRNADefWeight;

    // Default normalization parameters - Site level
    std::vector<float> mSiteDefLower;
    std::vector<float> mSiteDefUpper;
    std::vector<bool> mSiteDefReverse;

    // Default normalization parameters - RNA level
    std::vector<float> mRNADefLower;
    std::vector<float> mRNADefUpper;
    std::vector<bool> mRNADefReverse;

    std::string mTrue;
    std::string mFalse;

    std::string replace_site_key(std::string &pKey) {
        std::string key = pKey;  // copy
        key.replace(2, 1, "_");
        key.replace(7, 1, "_");

        return key;  // copy
    }

    std::string replace_rna_key(std::string &pKey) {
        std::string key = pKey;  // copy
        key.replace(2, 1, "_");
        key.replace(6, 1, "_");

        return key;  // copy
    }

    float get_float_from_map(const std::map<std::string, float> &pMap, std::string pKey) const {
        typedef std::map<std::string, float>::const_iterator TMFIt;
        TMFIt it = pMap.find(pKey);
        if (it != pMap.end()) {
            return (*it).second;
        } else {
            return false;
        }
    }

    bool get_bool_from_map(const std::map<std::string, bool> &pMap, std::string pKey) const {
        typedef std::map<std::string, bool>::const_iterator TMBIt;
        TMBIt it = pMap.find(pKey);
        if (it != pMap.end()) {
            return (*it).second;
        } else {
            return false;
        }
    }

};

} // namespace mkens

#endif //MKE_CONFIG_HPP_