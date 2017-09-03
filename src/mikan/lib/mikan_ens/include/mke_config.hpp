#ifndef MKE_CONFIG_HPP_
#define MKE_CONFIG_HPP_

#include <seqan/sequence.h>
#include "inih/INIReader.h"

namespace mkens {

struct ToolIdx {
    enum idx {
        MR, PT, RH, TM, TS, SV
    };
    static const unsigned Count = 6;
};

struct SiteIdx {
    enum idx {
        MRAlg, MREng,
        PTDdg, PTDpx, PTOpn,
        RHMfe,
        TSCtx,
        SVSvm
    };
    static const unsigned Count = 8;
};

struct RNAIdx {
    enum idx {
        MRAlg, MREng,
        PTDdg, PTDpx, PTOpn,
        RHMfe,
        TMSvm,
        TSCtx,
        SVSvm
    };
    static const unsigned Count = 9;
};

//
// Tool config interface
//
class MKEConfig {
public:
    static const int CommentPos = 26;

    // Define method
    MKEConfig() :
            mDESC("desc"),
            mASC("asc"),
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
        mRNAKeys[RNAIdx::PTDpx] = "pt:rna:dpx";
        mRNAKeys[RNAIdx::PTOpn] = "pt:rna:opn";
        mRNAKeys[RNAIdx::RHMfe] = "rh:rna:mfe";
        mRNAKeys[RNAIdx::TMSvm] = "tm:rna:svm";
        mRNAKeys[RNAIdx::TSCtx] = "ts:rna:ctx";
        mRNAKeys[RNAIdx::SVSvm] = "sv:rna:svm";

        // Score descriptions
        mRNADesc.resize(RNAIdx::Count);
        mRNADesc[RNAIdx::MRAlg] = "miRanda alignment score";
        mRNADesc[RNAIdx::MREng] = "miRanda energy score";
        mRNADesc[RNAIdx::PTDdg] = "PITA ddG";
        mRNADesc[RNAIdx::PTDpx] = "PITA dDuplex";
        mRNADesc[RNAIdx::PTOpn] = "PITA dOpen";
        mRNADesc[RNAIdx::RHMfe] = "RNAhybrid MFE";
        mRNADesc[RNAIdx::TMSvm] = "TargetMiner SVM";
        mRNADesc[RNAIdx::TSCtx] = "TargetScan context score";
        mRNADesc[RNAIdx::SVSvm] = "Two-step SVM score";

        // Score flags
        mRNADefFlg.resize(RNAIdx::Count);
        mRNADefFlg[RNAIdx::MRAlg] = true;
        mRNADefFlg[RNAIdx::MREng] = true;
        mRNADefFlg[RNAIdx::PTDdg] = true;
        mRNADefFlg[RNAIdx::PTDpx] = true;
        mRNADefFlg[RNAIdx::PTOpn] = true;
        mRNADefFlg[RNAIdx::RHMfe] = true;
        mRNADefFlg[RNAIdx::TMSvm] = true;
        mRNADefFlg[RNAIdx::TSCtx] = true;
        mRNADefFlg[RNAIdx::SVSvm] = true;

        //
        // Default weight values
        //
        // Site level
        mSiteDefWeight.resize(SiteIdx::Count);
        mSiteDefWeight[SiteIdx::MRAlg] = 1;
        mSiteDefWeight[SiteIdx::MREng] = 1;
        mSiteDefWeight[SiteIdx::PTDdg] = 1;
        mSiteDefWeight[SiteIdx::PTDpx] = 1;
        mSiteDefWeight[SiteIdx::PTOpn] = 1;
        mSiteDefWeight[SiteIdx::RHMfe] = 1;
        mSiteDefWeight[SiteIdx::TSCtx] = 1;
        mSiteDefWeight[SiteIdx::SVSvm] = 1;

        // RNA level
        mRNADefWeight.resize(RNAIdx::Count);
        mRNADefWeight[RNAIdx::MRAlg] = 1;
        mRNADefWeight[RNAIdx::MREng] = 1;
        mRNADefWeight[RNAIdx::PTDdg] = 1;
        mRNADefWeight[RNAIdx::PTDpx] = 1;
        mRNADefWeight[RNAIdx::PTOpn] = 1;
        mRNADefWeight[RNAIdx::RHMfe] = 1;
        mRNADefWeight[RNAIdx::TMSvm] = 1;
        mRNADefWeight[RNAIdx::TSCtx] = 1;
        mRNADefWeight[RNAIdx::SVSvm] = 1;

        //
        // Default normalization parameters - site level
        //
        mSiteDefLower.resize(SiteIdx::Count);
        mSiteDefUpper.resize(SiteIdx::Count);
        mSiteDefOrder.resize(SiteIdx::Count);

        // miRanda alignment score
        mSiteDefLower[SiteIdx::MRAlg] = 0;
        mSiteDefUpper[SiteIdx::MRAlg] = 1;
        mSiteDefOrder[SiteIdx::MRAlg] = mDESC;

        // miRanda energy score
        mSiteDefLower[SiteIdx::MREng] = 0;
        mSiteDefUpper[SiteIdx::MREng] = 1;
        mSiteDefOrder[SiteIdx::MREng] = mDESC;

        // PITA ddG
        mSiteDefLower[SiteIdx::PTDdg] = 0;
        mSiteDefUpper[SiteIdx::PTDdg] = 1;
        mSiteDefOrder[SiteIdx::PTDdg] = mDESC;

        // PITA dDuplex
        mSiteDefLower[SiteIdx::PTDpx] = 0;
        mSiteDefUpper[SiteIdx::PTDpx] = 1;
        mSiteDefOrder[SiteIdx::PTDpx] = mDESC;

        // PITA dOpen
        mSiteDefLower[SiteIdx::PTOpn] = 0;
        mSiteDefUpper[SiteIdx::PTOpn] = 1;
        mSiteDefOrder[SiteIdx::PTOpn] = mDESC;

        // RNAhybrid MFE
        mSiteDefLower[SiteIdx::RHMfe] = 0;
        mSiteDefUpper[SiteIdx::RHMfe] = 1;
        mSiteDefOrder[SiteIdx::RHMfe] = mDESC;

        // TargetScan context score
        mSiteDefLower[SiteIdx::TSCtx] = 0;
        mSiteDefUpper[SiteIdx::TSCtx] = 1;
        mSiteDefOrder[SiteIdx::TSCtx] = mDESC;

        // Two-step SVM score"
        mSiteDefLower[SiteIdx::SVSvm] = 0;
        mSiteDefUpper[SiteIdx::SVSvm] = 1;
        mSiteDefOrder[SiteIdx::SVSvm] = mDESC;

        //
        // Default normalization parameters - RNA level
        //
        // RNA score default parameters for normalization
        mRNADefLower.resize(RNAIdx::Count);
        mRNADefUpper.resize(RNAIdx::Count);
        mRNADefOrder.resize(RNAIdx::Count);

        // miRanda alignment score
        mRNADefLower[RNAIdx::MRAlg] = 0;
        mRNADefUpper[RNAIdx::MRAlg] = 1;
        mRNADefOrder[RNAIdx::MRAlg] = mDESC;

        // miRanda energy score
        mRNADefLower[RNAIdx::MREng] = 0;
        mRNADefUpper[RNAIdx::MREng] = 1;
        mRNADefOrder[RNAIdx::MREng] = mDESC;

        // PITA ddG
        mRNADefLower[RNAIdx::PTDdg] = 0;
        mRNADefUpper[RNAIdx::PTDdg] = 1;
        mRNADefOrder[RNAIdx::PTDdg] = mDESC;

        // PITA dDuplex
        mRNADefLower[RNAIdx::PTDpx] = 0;
        mRNADefUpper[RNAIdx::PTDpx] = 1;
        mRNADefOrder[RNAIdx::PTDpx] = mDESC;

        // PITA dOpen
        mRNADefLower[RNAIdx::PTOpn] = 0;
        mRNADefUpper[RNAIdx::PTOpn] = 1;
        mRNADefOrder[RNAIdx::PTOpn] = mDESC;

        // RNAhybrid MFE
        mRNADefLower[RNAIdx::RHMfe] = 0;
        mRNADefUpper[RNAIdx::RHMfe] = 1;
        mRNADefOrder[RNAIdx::RHMfe] = mDESC;

        // TargetMiner SVM
        mRNADefLower[RNAIdx::TMSvm] = 0;
        mRNADefUpper[RNAIdx::TMSvm] = 1;
        mRNADefOrder[RNAIdx::TMSvm] = mDESC;

        // TargetScan context score
        mRNADefLower[RNAIdx::TSCtx] = 0;
        mRNADefUpper[RNAIdx::TSCtx] = 1;
        mRNADefOrder[RNAIdx::TSCtx] = mDESC;

        // Two-step SVM score
        mRNADefLower[RNAIdx::SVSvm] = 0;
        mRNADefUpper[RNAIdx::SVSvm] = 1;
        mRNADefOrder[RNAIdx::SVSvm] = mDESC;

        init_config();
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

    float get_rna_lower(std::string &pKey) const {
        return get_float_from_map(mRNALower, pKey);
    }

    float get_rna_upper(std::string &pKey) const {
        return get_float_from_map(mRNAUpper, pKey);
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
    std::map<std::string, std::string> mSiteOrder;

    // Normalization parameters - RNA level
    std::map<std::string, float> mRNALower;
    std::map<std::string, float> mRNAUpper;
    std::map<std::string, std::string> mRNAOrder;

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
    std::vector<std::string> mSiteDefOrder;

    // Default normalization parameters - RNA level
    std::vector<float> mRNADefLower;
    std::vector<float> mRNADefUpper;
    std::vector<std::string> mRNADefOrder;

    std::string mDESC;
    std::string mASC;

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