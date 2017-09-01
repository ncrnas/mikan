#include <iostream>
#include <sstream>
#include "mke_config.hpp"        // MKEConfig

using namespace seqan;

namespace mkens {

//
// MKEConfig methods
//
void MKEConfig::init_config() {
    init_flags();
    init_weights();
    init_tool_config();
}

void MKEConfig::init_flags() {
    // Flags of tools
    for (unsigned i = 0; i < mToolKeys.size(); ++i) {
        mToolFlg[mToolKeys[i]] = mToolDefFlg[i];
    }

    // Flags of site level scores
    for (unsigned i = 0; i < mSiteKeys.size(); ++i) {
        mSiteFlg[mSiteKeys[i]] = mSiteDefFlg[i];
    }

    // Flags of RNA level scores
    for (unsigned i = 0; i < mRNAKeys.size(); ++i) {
        mRNAFlg[mRNAKeys[i]] = mRNADefFlg[i];
    }
}

void MKEConfig::init_weights() {
    // Weights for site level scores
    for (unsigned i = 0; i < mSiteKeys.size(); ++i) {
        mSiteWeight[mSiteKeys[i]] = mSiteDefWeight[i];
    }

    // Weights for RNA level scores
    for (unsigned i = 0; i < mRNAKeys.size(); ++i) {
        mRNAWeight[mRNAKeys[i]] = mRNADefWeight[i];
    }
}

void MKEConfig::init_tool_config() {
    // Site level config
    for (unsigned i = 0; i < mSiteKeys.size(); ++i) {
        mSiteLower[mSiteKeys[i]] = mSiteDefLower[i];
        mSiteUpper[mSiteKeys[i]] = mSiteDefUpper[i];
        mSiteOrder[mSiteKeys[i]] = mSiteDefOrder[i];
    }

    // RNA level config
    for (unsigned i = 0; i < mRNAKeys.size(); ++i) {
        mRNALower[mRNAKeys[i]] = mRNADefLower[i];
        mRNAUpper[mRNAKeys[i]] = mRNADefUpper[i];
        mRNAOrder[mRNAKeys[i]] = mRNADefOrder[i];
    }
}

void MKEConfig::print_config() {
    print_flags();
    print_weights();
    print_tool_config();
}

std::string &MKEConfig::get_bool_str(bool pFlag) {
    if (pFlag) {
        return mTrue;
    } else {
        return mFalse;
    }
}

void MKEConfig::print_flags() {
    std::cout << ";" << std::endl;
    std::cout << "; Flags that indicate if tools and scores be used in mikan" << std::endl;
    std::cout << ";" << std::endl;

    // Flags of tools
    std::cout << "[tool:flag]" << std::endl;
    for (unsigned i = 0; i < mToolKeys.size(); ++i) {
        std::cout << mToolKeys[i] << " = ";
        std::cout << std::left << std::setw(CommentPos - 6) << std::setfill(' ');
        std::cout << get_bool_str(mToolFlg[mToolKeys[i]]);
        std::cout << "; " << mToolDesc[i] << std::endl;
    }
    std::cout << std::endl;

    // Flags of site level scores
    std::cout << "[site:flag]" << std::endl;
    for (unsigned i = 0; i < mSiteKeys.size(); ++i) {
        std::string ckey = replace_site_key(mSiteKeys[i]);
        std::cout << ckey << " = ";
        std::cout << std::left << std::setw(CommentPos - 15) << std::setfill(' ');
        std::cout << get_bool_str(mSiteFlg[mSiteKeys[i]]);
        std::cout << "; " << mSiteDesc[i] << std::endl;
    }
    std::cout << std::endl;

    // Flags of RNA level scores
    std::cout << "[rna:flag]" << std::endl;
    for (unsigned i = 0; i < mRNAKeys.size(); ++i) {
        std::string ckey = replace_rna_key(mRNAKeys[i]);
        std::cout << ckey << " = ";
        std::cout << std::left << std::setw(CommentPos - 14) << std::setfill(' ');
        std::cout << get_bool_str(mRNAFlg[mRNAKeys[i]]);
        std::cout << "; " << mRNADesc[i] << std::endl;
    }
    std::cout << std::endl;
}

void MKEConfig::print_weights() {
    std::cout << ";" << std::endl;
    std::cout << "; Weight values to calculate final mikan score" << std::endl;
    std::cout << ";" << std::endl;

    // Weights for site level scores
    std::cout << "[site:weight]" << std::endl;
    for (unsigned i = 0; i < mSiteKeys.size(); ++i) {
        std::string ckey = replace_site_key(mSiteKeys[i]);
        std::cout << ckey << " = ";
        std::cout << std::left << std::setw(CommentPos - 15) << std::setfill(' ');
        std::cout << mSiteWeight[mSiteKeys[i]];
        std::cout << "; " << mSiteDesc[i] << std::endl;
    }
    std::cout << std::endl;

    // Weights for RNA level scores
    std::cout << "[rna:weight]" << std::endl;
    for (unsigned i = 0; i < mRNAKeys.size(); ++i) {
        std::string ckey = replace_rna_key(mRNAKeys[i]);
        std::cout << ckey << " = ";
        std::cout << std::left << std::setw(CommentPos - 14) << std::setfill(' ');
        std::cout << mRNAWeight[mRNAKeys[i]];
        std::cout << "; " << mRNADesc[i] << std::endl;
    }
    std::cout << std::endl;
}

void MKEConfig::print_tool_config() {
    std::cout << ";" << std::endl;
    std::cout << "; Parameters used to normalize scores to range [0, 1] - site level" << std::endl;
    std::cout << ";" << std::endl;

    // Site level config
    for (unsigned i = 0; i < mSiteKeys.size(); ++i) {
        std::cout << "; " << mSiteDesc[i] << std::endl;

        std::cout << "[" << mSiteKeys[i] << "]" << std::endl;

        std::cout << "lower = ";
        std::cout << std::left << std::setw(CommentPos - 9) << std::setfill(' ');
        std::cout << mSiteLower[mSiteKeys[i]];
        std::cout << "; Lower bound" << std::endl;

        std::cout << "upper = ";
        std::cout << std::left << std::setw(CommentPos - 9) << std::setfill(' ');
        std::cout << mSiteUpper[mSiteKeys[i]];
        std::cout << "; Upper bound" << std::endl;

        std::cout << "order = ";
        std::cout << std::left << std::setw(CommentPos - 9) << std::setfill(' ');
        std::cout << mSiteOrder[mSiteKeys[i]];
        std::cout << "; Score order (desc or asc)" << std::endl;

        std::cout << std::endl;
    }

    std::cout << ";" << std::endl;
    std::cout << "; Parameters used to normalize scores to range [0, 1] - RNA level" << std::endl;
    std::cout << ";" << std::endl;

    // RNA level config
    for (unsigned i = 0; i < mRNAKeys.size(); ++i) {
        std::cout << "; " << mRNADesc[i] << std::endl;

        std::cout << "[" << mRNAKeys[i] << "]" << std::endl;

        std::cout << "lower = ";
        std::cout << std::left << std::setw(CommentPos - 9) << std::setfill(' ');
        std::cout << mRNALower[mRNAKeys[i]];
        std::cout << "; Lower bound" << std::endl;

        std::cout << "upper = ";
        std::cout << std::left << std::setw(CommentPos - 9) << std::setfill(' ');
        std::cout << mRNAUpper[mRNAKeys[i]];
        std::cout << "; Upper bound" << std::endl;

        std::cout << "order = ";
        std::cout << std::left << std::setw(CommentPos - 9) << std::setfill(' ');
        std::cout << mRNAOrder[mRNAKeys[i]];
        std::cout << "; Score order (desc or asc)" << std::endl;

        std::cout << std::endl;
    }
}

int MKEConfig::parse_config(std::string &pConfFile) {

    INIReader reader(pConfFile);
    if (reader.ParseError() < 0) {
        return 1;
    }

    int retVal;
    retVal = parse_flags(reader);
    if (retVal != 0) {
        return 1;
    }

    retVal = parse_weights(reader);
    if (retVal != 0) {
        return 1;
    }

    retVal = parse_tool_config(reader);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

int MKEConfig::parse_flags(INIReader &pReader) {
    // Flags of tools
    for (TKeyIt it = mToolKeys.begin(); it != mToolKeys.end(); ++it) {
        mToolFlg[*it] = pReader.GetBoolean("tool:flag", *it, mToolFlg[*it]);
    }

    // Flags of site level scores
    for (TKeyIt it = mSiteKeys.begin(); it != mSiteKeys.end(); ++it) {
        std::string ckey = replace_site_key(*it);
        mSiteFlg[*it] = pReader.GetBoolean("site:flag", ckey, mSiteFlg[*it]);
    }

    // Flags of RNA level scores
    for (TKeyIt it = mRNAKeys.begin(); it != mRNAKeys.end(); ++it) {
        std::string ckey = replace_rna_key(*it);
        mRNAFlg[*it] = pReader.GetBoolean("rna:flag", ckey, mRNAFlg[*it]);
    }

    return 0;
}

int MKEConfig::parse_weights(INIReader &pReader) {
    // Weights for site level scores
    for (TKeyIt it = mSiteKeys.begin(); it != mSiteKeys.end(); ++it) {
        std::string ckey = replace_site_key(*it);
        mSiteWeight[*it] = static_cast<float>(pReader.GetReal("site:weight", ckey, mSiteWeight[*it]));
    }

    // Weights for RNA level scores
    for (TKeyIt it = mRNAKeys.begin(); it != mRNAKeys.end(); ++it) {
        std::string ckey = replace_rna_key(*it);
        mRNAWeight[*it] = static_cast<float>(pReader.GetReal("rna:weight", ckey, mRNAWeight[*it]));
    }

    return 0;
}

int MKEConfig::parse_tool_config(INIReader &pReader) {

    // Site level config
    for (TKeyIt it = mSiteKeys.begin(); it != mSiteKeys.end(); ++it) {
        mSiteLower[*it] = static_cast<float>(pReader.GetReal(*it, "lower", mSiteLower[*it]));

        mSiteUpper[*it] = static_cast<float>(pReader.GetReal(*it, "upper", mSiteUpper[*it]));

        mSiteOrder[*it] = pReader.Get(*it, "order", mSiteOrder[*it]);
    }

    // RNA level config
    for (TKeyIt it = mRNAKeys.begin(); it != mRNAKeys.end(); ++it) {
        mRNALower[*it] = static_cast<float>(pReader.GetReal(*it, "lower", mRNALower[*it]));

        mRNAUpper[*it] = static_cast<float>(pReader.GetReal(*it, "upper", mRNAUpper[*it]));

        mRNAOrder[*it] = pReader.Get(*it, "order", mRNAOrder[*it]);
    }

    return 0;
}

} // namespace mkens

