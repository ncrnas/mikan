#ifndef VR16_DDG_CORE_HPP_
#define VR16_DDG_CORE_HPP_

#include <vr16_duplex.hpp>                // VR16Duplex
#include <vr16_energy.hpp>                // VR16EnergyParams
#include <vr16_fold.hpp>                  // VR16Fold
#include <vr16_fold_options.hpp>          // VR16FoldOptions
#include <vr16_pair_mat.hpp>              // VR16PairMat
#include <vr16_params.hpp>                // VR16Params, VR16ParamIL, VR16PFParams
#include <vr16_part_func.hpp>             // VR16PartFunc
#include <vector>
#include <string>

namespace vr16 {

//
// Return values for ViennaRNA duplexfold
//
class VR16DuplexRet {
public:
    // Declare variables
    std::vector<std::string> mStructure;
    std::vector<int> mL1;
    std::vector<int> mStartMiR;
    std::vector<int> mEndMiR;
    std::vector<int> mStartTarget;
    std::vector<int> mEndTarget;
    std::vector<int> mMiRLen;
    std::vector<int> mTargetLen;
    std::vector<double> mDGall;
    std::vector<double> mDG5;
    std::vector<double> mDG3;

public:
    // Define methods
    VR16DuplexRet() : mVecSize(0) {}

    // Method prototypes
    void init_ret_vals(int pSize);

private:
    int mVecSize;

};

//
// Return values for ddG4 calculation
//
class VR16DDG4Ret {
public:
    // Constant values
    static const int ARRSIZE = 50;

    // Declare variables
    std::vector<double> mDG0;
    std::vector<std::vector<double> > mD1Array;
    std::vector<std::vector<double> > mD2Array;
    std::vector<std::vector<double> > mDDGArray;
    std::vector<double> mDDGSum;
    std::vector<double> mP;
    int mNExps;

public:
    // Define methods
    VR16DDG4Ret() : mNExps(0), mVecSize(0) {}

    // Method prototypes
    void init_ret_vals(int pSize);

private:
    int mVecSize;

};

//
// Main ddG calculation
//
class VR16DDGWorkSpace {
public:
    // Constant values
    const double MIN_EXP_DIFF;

    // Declare variables
    double mTemperature;

public:
    VR16DDGWorkSpace() :
            MIN_EXP_DIFF(-100.0), mTemperature(37.0), SFACT(1.07), mEn(), mParams(mEn, mTemperature),
            mPairMat(), mDupOpts(), mDG4Opts(), mPFParams(mEn, mTemperature, mParams, mDG4Opts),
            mPpIL(mEn, mParams, mPairMat, mDG4Opts),
            mDup(mEn, mParams, mPpIL, mPairMat, mDupOpts), mFold(mEn, mParams, mPpIL, mPairMat, mDG4Opts),
            mPf(mEn, mParams, mPFParams, mPairMat, mDG4Opts), mDupRet(), mDDG4Ret() {
        init_workspace();
    }

    void set_duplex_backtrack(bool pBT) { mDupOpts.mDoBacktrack = pBT; }

    void preppare_duplexfold(int pSize);

    int duplexfold(int pRetIdx, std::string &pS1, std::string &pS2, std::vector<int> &pArrayI,
                   std::vector<int> &pArrayJ);

    void print_duplexfold_ret_vals(int pRetIdx);

    void prepare_ddg4(int pSize, int pUpRest, int pTargetStart, int pRestrictedFrom, int pRestrictedTo);

    int calc_ddg4(int pRetIdx, std::string &pString);

    void print_ddg4_ret_vals(int pRetIdx);

    void print_array_size(int pI);

    double get_dgall(int pRetIdx) { return mDupRet.mDGall[pRetIdx]; }

    double get_dg5(int pRetIdx) { return mDupRet.mDG5[pRetIdx]; }

    double get_dg3(int pRetIdx) { return mDupRet.mDG3[pRetIdx]; }

    double get_dg0(int pRetIdx) { return mDDG4Ret.mDG0[pRetIdx]; }

    double get_dg1(int pRetIdx) { return mDDG4Ret.mD1Array[pRetIdx][0]; }

    double get_dgsum(int pRetIdx) { return mDDG4Ret.mDDGSum[pRetIdx]; }

    const std::string &get_structure(int pRetIdx) { return mDupRet.mStructure[pRetIdx]; }

    int get_l1(int pRetIdx) { return mDupRet.mL1[pRetIdx]; }

private:
    void init_workspace();

    double log_of_sum_of_exps(std::vector<double> &pExps, int pNumExps);

private:
    const double SFACT;

    VR16EnergyParams mEn;
    VR16Params mParams;
    VR16PairMat mPairMat;
    VR16FoldOptions mDupOpts;
    VR16FoldOptions mDG4Opts;
    VR16PFParams mPFParams;
    VR16ParamIL mPpIL;

    VR16Duplex mDup;
    VR16Fold mFold;
    VR16PartFunc mPf;

    VR16DuplexRet mDupRet;
    VR16DDG4Ret mDDG4Ret;

    int mUpRest;
    int mTargetStart;
    int mRestrictedFrom;
    int mRestrictedTo;

};

} // namespace vr16

#endif /* VR16_DDG_CORE_HPP_ */
